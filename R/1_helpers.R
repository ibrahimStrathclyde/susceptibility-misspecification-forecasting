# =============================================================================
# 1_helpers.R
#
# Shared functions for the susceptibility-distribution misspecification analyses.
# The file contains discretisation functions, SEIR ODE systems, simulation
# routines, Poisson likelihoods, and maximum-likelihood fitting helpers used by
# the analysis scripts in this repository.
# =============================================================================

library(deSolve)
library(dplyr)

# =============================================================================
# 1.  UTILITY FUNCTIONS
# =============================================================================

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))

# =============================================================================
# 2.  INITIAL-CONDITION BUILDER
# =============================================================================

initDistr <- function(q, K, E0, I0, N) {
  if (trunc(K / 2) * 2 != K) stop("K must be even")

  proportion_E <- E0 / (E0 + I0)
  epsilon      <- E0 / (N * proportion_E)

  init_s <- (1 - epsilon) * q * N
  init_e <- as.vector(rmultinom(1, E0, q))
  init_i <- as.vector(rmultinom(1, I0, q))
  init_r <- rep(0, K)
  init_c <- rep(0, K)

  names(init_s) <- paste0("S", 1:K)
  names(init_e) <- paste0("E", 1:K)
  names(init_i) <- paste0("I", 1:K)
  names(init_r) <- paste0("R", 1:K)
  names(init_c) <- paste0("C", 1:K)

  init <- c(round(init_s), init_e, init_i, init_r, init_c)
  return(init)
}

# =============================================================================
# 3.  HETEROGENEOUS SEIR ODE
# =============================================================================

hetsus_model.ct <- function(t, y, parms) {
  x       <- parms[grepl("^x", names(parms))]
  parms_2 <- parms[!grepl("^x", names(parms))]
  z_list  <- as.list(parms_2)
  z_list$x <- x;  z_list$t <- t;  z_list$y <- y

  with(z_list, {
    S <- y[1:K]
    E <- y[(K + 1):(2 * K)]
    I <- y[(2 * K + 1):(3 * K)]
    R <- y[(3 * K + 1):(4 * K)]
    C <- y[(4 * K + 1):(5 * K)]

    # Piecewise-linear NPI profile
    if (t <= t0) {
      prox <- c_value1
    } else if (t <= t1) {
      prox <- c_value1 - (c_value1 - c_value2) * (t - t0) / (t1 - t0)
    } else if (t <= t2) {
      prox <- c_value2
    } else if (t <= t3) {
      prox <- c_value2 + (c_value3 - c_value2) * (t - t2) / (t3 - t2)
    } else {
      prox <- c_value3
    }

    Beta   <- R0 * prox / (rho / delta + 1 / gamma)
    lambda <- Beta / N * (rho * sum(E) + sum(I))

    dS <- -lambda * x * S
    dE <-  lambda * x * S - delta * E
    dI <-  delta * E - gamma * I
    dR <-  gamma * I
    dC <-  delta * E

    return(list(c(dS, dE, dI, dR, dC)))
  })
}

# =============================================================================
# 4.  DISCRETISATION FUNCTIONS  
# =============================================================================

# ---------- 4a.  GAMMA ----------
Discretize_gamma_LA_SM <- function(
    n_groups,
    alpha,
    beta         = alpha,
    spacing      = c("equal", "linear"),
    rep          = c("condexp", "midpoint"),
    calibration  = c("log_affine", "split_mult"),
    positivity_guard       = TRUE,
    min_positive_threshold = 1e-8,
    print_message          = FALSE
) {
  if (alpha <= 0) stop("Alpha must be positive")
  if (beta  <= 0) stop("Beta must be positive")

  spacing     <- match.arg(spacing)
  rep         <- match.arg(rep)
  calibration <- match.arg(calibration)

  m          <- alpha / beta          # target mean
  target_var <- alpha / (beta^2)      # target variance

  # --- spacing ---
  if (spacing == "linear") {
    n_g <- trunc(n_groups / 2)
    q1 <- 1:n_g;  q2 <- n_g:1
    q <- c(q1, q2);  q <- q / sum(q)
  } else {
    q <- rep(1 / n_groups, n_groups)
  }

  # --- quantile edges & representatives ---
  z      <- c(0, cumsum(q))
  edges  <- qgamma(z, shape = alpha, rate = beta)
  a      <- edges[-length(edges)]
  b      <- edges[-1]
  mid_p  <- (z[-1] + z[-(n_groups + 1)]) / 2
  x_mid  <- qgamma(mid_p, shape = alpha, rate = beta)

  P_bin  <- pgamma(b, shape = alpha, rate = beta) -
            pgamma(a, shape = alpha, rate = beta)
  E1_bin <- (alpha / beta) *
            (pgamma(b, shape = alpha + 1, rate = beta) -
             pgamma(a, shape = alpha + 1, rate = beta))
  x_condexp <- E1_bin / pmax(P_bin, .Machine$double.xmin)

  if (rep == "midpoint") x <- x_mid else x <- x_condexp
  x_init <- x

  # --- mean correction on last bin ---
  xbar <- sum(q * x)
  x[n_groups] <- (m - xbar + x[n_groups] * q[n_groups]) / q[n_groups]

  # --- solution_text for backward compatibility with objective-function checks ---
  solution_text      <- "log_affine"
  calibration_params <- NULL

  # ====== LOG-AFFINE CALIBRATION ======
  if (calibration == "log_affine") {
    logsumexp_w <- function(v) { M <- max(v); M + log(sum(exp(v - M))) }
    h_of_t      <- function(t, y, q) logsumexp_w(log(q) + t * y)
    y <- log(if (positivity_guard) pmax(x, min_positive_threshold) else x)

    rhs  <- log(1 + 1 / alpha)
    g_B  <- function(B) h_of_t(2 * B, y, q) - 2 * h_of_t(B, y, q) - rhs

    Blo <- 0;  Bhi <- 1;  ghi <- g_B(Bhi);  it <- 0
    while (ghi <= 0 && Bhi < 100 && it < 50) {
      Bhi <- 2 * Bhi + 1e-9;  ghi <- g_B(Bhi);  it <- it + 1
    }

    if (ghi <= 0) {
      solution_text      <- "log_affine (no bracket)"
      calibration_params <- list(type = "log_affine_failed")
    } else {
      B <- uniroot(g_B, interval = c(Blo, Bhi), tol = 1e-12)$root
      A <- log(m) - h_of_t(B, y, q)
      x <- exp(A + B * y)
      x <- x / sum(q * x) * m
      solution_text      <- "log_affine"
      calibration_params <- list(type = "log_affine", A = A, B = B,
                                 mean_check = sum(q * x),
                                 var_check  = sum(q * (x - sum(q * x))^2))
    }

  # ====== SPLIT-MULTIPLICATIVE CALIBRATION ======
  } else if (calibration == "split_mult") {
    if (trunc(n_groups / 2) * 2 != n_groups) stop("split_mult requires even n_groups")
    ng2 <- n_groups / 2
    L <- 1:ng2;  R <- (ng2 + 1):n_groups
    X1 <- sum(q[L] * x[L]);  X2 <- sum(q[R] * x[R])
    M1 <- sum(q[L] * x[L]^2);  M2 <- sum(q[R] * x[R]^2)
    Treq <- m^2 + target_var

    a_coef <- M1 + M2 * (X1^2) / (X2^2)
    b_coef <- -2 * M2 * m * X1 / (X2^2)
    c_coef <- M2 * (m^2) / (X2^2) - Treq
    disc   <- b_coef^2 - 4 * a_coef * c_coef

    multiplier_cap  <- 3
    fallback_needed <- (!is.finite(disc) || disc <= 0 ||
                        !is.finite(a_coef) || a_coef == 0)
    mu1_best <- mu2_best <- NA_real_;  x_best <- NULL

    if (!fallback_needed) {
      roots <- c((-b_coef + sqrt(disc)) / (2 * a_coef),
                 (-b_coef - sqrt(disc)) / (2 * a_coef))
      best_score <- Inf
      for (yc in roots) {
        mu1 <- 1 / yc;  mu2 <- (m - X1 / mu1) / X2
        if (!is.finite(mu1) || !is.finite(mu2) || mu1 <= 0 || mu2 <= 0) next
        if (abs(log(mu1)) > multiplier_cap || abs(log(mu2)) > multiplier_cap) next
        cand <- x;  cand[L] <- cand[L] / mu1;  cand[R] <- cand[R] * mu2
        if (max(cand[L]) > min(cand[R])) next
        sc <- abs(log(mu1)) + abs(log(mu2))
        if (sc < best_score) {
          best_score <- sc;  mu1_best <- mu1;  mu2_best <- mu2;  x_best <- cand
        }
      }
      if (!is.finite(mu1_best) || !is.finite(mu2_best)) fallback_needed <- TRUE
    }

    if (fallback_needed) {
      if (print_message) cat("split_mult: no admissible root; identity multipliers.\n")
      solution_text      <- "split_mult_identity"
      calibration_params <- list(type = "split_mult_identity", mu1 = 1, mu2 = 1)
    } else {
      x <- x_best
      var_now <- sum(q * (x - m)^2)
      if (is.finite(var_now) && var_now > 0 && abs(var_now - target_var) > 1e-10) {
        s <- sqrt(target_var / var_now);  x <- m + s * (x - m);  x <- x / sum(q * x) * m
      } else {
        x <- x / sum(q * x) * m
      }
      solution_text      <- "split_mult"
      calibration_params <- list(type = "split_mult", mu1 = mu1_best, mu2 = mu2_best)
    }
  }

  # --- summary ---
  mean_emp <- sum(q * x);  var_emp <- sum(q * (x - mean_emp)^2)
  cv_emp   <- sqrt(var_emp) / mean_emp

  list(
    q = q,  x = x,  x_init = x_init,  x_mid = x_mid,
    edges = edges,  a = a,  b = b,
    solution_text = solution_text,
    calibration_params = calibration_params,
    summary = list(
      empirical   = c(mean = mean_emp, variance = var_emp, CV = cv_emp),
      closed_form = c(mean = m, variance = target_var, CV = 1 / sqrt(alpha))
    )
  )
}


# ---------- 4b.  LOGNORMAL ----------
Discretize_lognormal_LA_SM <- function(
    n_groups,
    CV,
    spacing      = c("equal", "linear"),
    rep          = c("condexp", "midpoint"),
    calibration  = c("log_affine", "split_mult"),
    positivity_guard       = TRUE,
    min_positive_threshold = 1e-3,
    print_message          = FALSE
) {
  spacing     <- match.arg(spacing)
  rep         <- match.arg(rep)
  calibration <- match.arg(calibration)
  stopifnot(n_groups >= 1, CV > 0)

  sigma      <- sqrt(log(1 + CV^2))
  mu         <- -0.5 * sigma^2
  var_target <- CV^2

  # --- spacing ---
  if (spacing == "equal") {
    q <- rep(1 / n_groups, n_groups)
  } else {
    half_n <- floor(n_groups / 2)
    q1 <- 1:half_n;  q2 <- rev(q1)
    if (n_groups %% 2 == 1) q <- c(q1, max(q1), q2) else q <- c(q1, q2)
    q <- q / sum(q)
  }

  z     <- c(0, cumsum(q))
  edges <- qlnorm(z, meanlog = mu, sdlog = sigma)
  a     <- edges[-length(edges)]
  b     <- edges[-1]
  mid_p <- (z[-1] + z[-length(z)]) / 2
  x_mid <- qlnorm(mid_p, meanlog = mu, sdlog = sigma)

  if (rep == "midpoint") {
    x_rep <- x_mid
  } else {
    P_bin  <- plnorm(b, meanlog = mu, sdlog = sigma) -
              plnorm(a, meanlog = mu, sdlog = sigma)
    u_b    <- (log(b) - mu - sigma^2) / sigma
    u_a    <- (log(a) - mu - sigma^2) / sigma
    E1_bin <- exp(mu + 0.5 * sigma^2) * (pnorm(u_b) - pnorm(u_a))
    x_rep  <- E1_bin / pmax(P_bin, .Machine$double.xmin)
  }

  x0     <- x_rep / sum(q * x_rep)   # mean-normalise
  x_init <- x0;  x <- x0

  solution_text      <- "log_affine"
  calibration_params <- NULL

  # ====== LOG-AFFINE CALIBRATION ======
  if (calibration == "log_affine") {
    logsumexp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
    h_of_t    <- function(t, y, q) logsumexp(log(q) + t * y)
    y <- log(if (positivity_guard) pmax(x, min_positive_threshold) else x)

    g_of_B <- function(B) h_of_t(2 * B, y, q) - 2 * h_of_t(B, y, q) - log(1 + CV^2)

    B_lo <- 0;  B_hi <- 1;  g_hi <- g_of_B(B_hi);  iter <- 0
    while (g_hi <= 0 && B_hi < 100 && iter < 50) {
      B_hi <- 2 * B_hi + 1e-9;  g_hi <- g_of_B(B_hi);  iter <- iter + 1
    }

    if (g_hi <= 0) {
      if (print_message) cat("log_affine: could not bracket root; keeping initial.\n")
      solution_text      <- "log_affine (no bracket)"
      calibration_params <- list(type = "log_affine_identity")
      x <- x_init
    } else {
      B_sol <- uniroot(g_of_B, interval = c(B_lo, B_hi), tol = 1e-12)$root
      A_sol <- -h_of_t(B_sol, y, q)
      x     <- exp(A_sol + B_sol * y)
      x     <- x / sum(q * x)
      solution_text      <- "log_affine"
      calibration_params <- list(type = "log_affine", A = A_sol, B = B_sol)
    }

  # ====== SPLIT-MULTIPLICATIVE CALIBRATION ======
  } else {
    if (n_groups %% 2L != 0L) stop("split_mult requires even n_groups")
    ng2 <- n_groups / 2
    L <- 1:ng2;  R <- (ng2 + 1):n_groups
    X1 <- sum(q[L] * x[L]);  X2 <- sum(q[R] * x[R])
    M1 <- sum(q[L] * x[L]^2);  M2 <- sum(q[R] * x[R]^2)
    Treq <- 1 + var_target

    a_coef <- M1 + M2 * (X1^2) / (X2^2)
    b_coef <- -2 * M2 * X1 / (X2^2)
    c_coef <- M2 / (X2^2) - Treq
    disc   <- b_coef^2 - 4 * a_coef * c_coef

    multiplier_cap  <- 3
    fallback_needed <- (!is.finite(disc) || disc <= 0 ||
                        !is.finite(a_coef) || a_coef == 0)
    mu1 <- mu2 <- NA_real_;  x_sm <- NULL

    if (!fallback_needed) {
      roots <- c((-b_coef + sqrt(disc)) / (2 * a_coef),
                 (-b_coef - sqrt(disc)) / (2 * a_coef))
      best_score <- Inf
      for (r in roots) {
        mu1_try <- r;  mu2_try <- (1 - mu1_try * X1) / X2
        if (!is.finite(mu1_try) || !is.finite(mu2_try) ||
            mu1_try <= 0 || mu2_try <= 0) next
        if (abs(log(mu1_try)) > multiplier_cap ||
            abs(log(mu2_try)) > multiplier_cap) next
        cand <- x;  cand[L] <- mu1_try * cand[L];  cand[R] <- mu2_try * cand[R]
        if (max(cand[L]) > min(cand[R])) next
        cand <- cand / sum(q * cand)
        sc <- abs(log(mu1_try)) + abs(log(mu2_try))
        if (sc < best_score) {
          best_score <- sc;  mu1 <- mu1_try;  mu2 <- mu2_try;  x_sm <- cand
        }
      }
      if (!is.finite(mu1) || !is.finite(mu2) || is.null(x_sm)) fallback_needed <- TRUE
    }

    if (fallback_needed) {
      if (print_message) cat("split_mult: no admissible root; identity multipliers.\n")
      solution_text      <- "split_mult_identity"
      calibration_params <- list(type = "split_mult_identity", mu1 = 1, mu2 = 1)
      x <- x_init
    } else {
      x <- x_sm
      solution_text      <- "split_mult"
      calibration_params <- list(type = "split_mult", mu1 = mu1, mu2 = mu2)
    }
  }

  mean_emp <- sum(q * x);  var_emp <- sum(q * (x - mean_emp)^2)
  cv_emp   <- sqrt(var_emp) / mean_emp

  list(
    q = q,  x = x,  x_init = x_init,  x_mid = x_mid,
    edges = edges,  a = a,  b = b,
    solution_text = solution_text,
    calibration_params = calibration_params,
    summary = list(
      empirical   = c(mean = mean_emp, variance = var_emp, CV = cv_emp),
      closed_form = c(mean = 1, variance = var_target, CV = CV)
    )
  )
}


# =============================================================================
# 5.  SIMULATION FUNCTIONS
# =============================================================================

simulate_cases_hetsus_model_gamma <- function(
    alpha, K, R0, delta, rho, gamma, N, E0, I0,
    t0, t1, t2, t3, c_value1 = 1, c_value2, c_value3 = 1, tfinal
) {
  xval <- Discretize_gamma_LA_SM(
    n_groups = K, alpha = alpha, beta = alpha,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )

  init_het <- initDistr(q = xval$q, K = K, E0 = E0, I0 = I0, N = N)

  parms <- c(xval$x, K, R0, gamma, rho, delta, N,
             c_value1, c_value2, c_value3, t0, t1, t2, t3)
  names(parms) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  times <- if (exists("times_full")) times_full else 0:tfinal

  outdf <- as.data.frame(ode(y = init_het, times = times,
                              func = hetsus_model.ct, parms = parms))

  C   <- rowSums(outdf[, (4 * K + 2):(5 * K + 1)], na.rm = TRUE)
  inc <- c(0, diff(C))
  lam <- pmax(inc, 0)

  cases    <- rpois(length(lam), lambda = lam)
  sim_data <- data.frame(time = times[-1], reports = cases[-1])
  return(list(sim_data = sim_data))
}


simulate_cases_hetsus_model_lognormal <- function(
    CV, K, R0, delta, rho, gamma, N, E0, I0,
    t0, t1, t2, t3, c_value1 = 1, c_value2, c_value3 = 1, tfinal
) {
  xval <- Discretize_lognormal_LA_SM(
    n_groups = K, CV = CV,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )

  init_het <- initDistr(q = xval$q, K = K, E0 = E0, I0 = I0, N = N)

  parms <- c(xval$x, K, R0, gamma, rho, delta, N,
             c_value1, c_value2, c_value3, t0, t1, t2, t3)
  names(parms) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  times <- if (exists("times_full")) times_full else 0:tfinal

  outdf <- as.data.frame(ode(y = init_het, times = times,
                              func = hetsus_model.ct, parms = parms))

  C   <- rowSums(outdf[, (4 * K + 2):(5 * K + 1)], na.rm = TRUE)
  inc <- c(0, diff(C))
  lam <- pmax(inc, 0)

  cases    <- rpois(length(lam), lambda = lam)
  sim_data <- data.frame(time = times[-1], reports = cases[-1])
  return(list(sim_data = sim_data))
}


# =============================================================================
# 6.  POISSON LOG-LIKELIHOOD  (gamma and lognormal paths — identical logic)
# =============================================================================

poisson.loglik <- function(params, sim.data, init_het, times) {
  outdf <- as.data.frame(ode(y = init_het, times = times_full,
                              func = hetsus_model.ct, parms = params))
  C   <- rowSums(outdf[, (4 * K + 2):(5 * K + 1)], na.rm = TRUE)
  inc <- c(0, diff(C))
  df  <- data.frame(time = outdf[, 1], Inc = inc)

  if (!is.data.frame(sim.data)) sim.data <- as.data.frame(sim.data)
  obs_col <- if ("reports" %in% names(sim.data)) "reports" else "cases"
  if (!obs_col %in% names(sim.data)) stop("No 'reports' or 'cases' column.")

  df <- df %>% filter(time %in% sim.data$time)
  if (nrow(df) != nrow(sim.data)) stop("Time mismatch between model and data.")

  lambda_ <- ifelse(df$Inc == 0, 0.0001, df$Inc)
  sum(dpois(x = sim.data[[obs_col]], lambda = lambda_, log = TRUE))
}

# Lognormal path — identical logic, kept separate for clarity
poisson.loglik.lognormal <- function(params, sim.data, init_het, times) {
  poisson.loglik(params, sim.data, init_het, times)
}


# =============================================================================
# 7.  SINGLE-EPIDEMIC OBJECTIVE FUNCTIONS
# =============================================================================

# --- Gamma ---
f4_optim_1epic_poisson.loglikwithNPI <- function(par, sim.data, init_het) {
  z_R0 <- exp(par[1]);  z_v <- exp(par[2])
  z_t0 <- exp(par[3]);  z_c2 <- expit(par[4])
  if (!all(is.finite(c(z_R0, z_v, z_t0, z_c2)))) return(.Machine$double.xmax)

  alpha_fit <- 1 / (z_v^2)
  xval <- Discretize_gamma_LA_SM(
    n_groups = K, alpha = as.numeric(alpha_fit), beta = as.numeric(alpha_fit),
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  if (xval$solution_text != "log_affine") return(.Machine$double.xmax)
  if (is.null(xval$x) || length(xval$x) != K ||
      any(!is.finite(xval$x)) || any(xval$x <= 0)) return(.Machine$double.xmax)

  params <- c(xval$x,
              K = K, R0 = z_R0, gamma = gamma_spec, rho = rho_spec,
              delta = delta_spec, N = N,
              c_value1 = c_value1_spec, c_value2 = z_c2, c_value3 = c_value3_spec,
              t0 = z_t0, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec)
  names(params) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  tryCatch(
    -poisson.loglik(params, sim.data = sim.data, init_het = init_het, times = sim.data$time),
    error = function(e) .Machine$double.xmax
  )
}

# --- Lognormal ---
f4_optim_1epic_poisson.loglikwithNPI_lognormal <- function(par, sim.data, init_het) {
  z_R0 <- exp(par[1]);  z_cv <- exp(par[2])
  z_t0 <- exp(par[3]);  z_c2 <- expit(par[4])
  if (!all(is.finite(c(z_R0, z_cv, z_t0, z_c2)))) return(.Machine$double.xmax)

  xval <- Discretize_lognormal_LA_SM(
    n_groups = K, CV = z_cv,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  if (xval$solution_text != "log_affine") return(.Machine$double.xmax)
  if (is.null(xval$x) || length(xval$x) != K ||
      any(!is.finite(xval$x)) || any(xval$x <= 0)) return(.Machine$double.xmax)

  params <- c(xval$x,
              K = K, R0 = z_R0, gamma = gamma_spec, rho = rho_spec,
              delta = delta_spec, N = N,
              c_value1 = c_value1_spec, c_value2 = z_c2, c_value3 = c_value3_spec,
              t0 = z_t0, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec)
  names(params) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  tryCatch(
    -poisson.loglik.lognormal(params, sim.data = sim.data, init_het = init_het, times = sim.data$time),
    error = function(e) .Machine$double.xmax
  )
}


# =============================================================================
# 8.  SINGLE-EPIDEMIC FITTING WRAPPERS
# =============================================================================

fit4_seir_1epic_poisson.loglikwithNPI <- function(dat, init_het) {
  start_par <- c(log(2), log(1.1), log(12), logit(0.2))

  fit1 <- optim(
    par = start_par,
    fn  = f4_optim_1epic_poisson.loglikwithNPI,
    sim.data = dat, init_het = init_het,
    method  = "Nelder-Mead",
    control = list(trace = 0, maxit = 800),
    hessian = TRUE
  )

  H <- fit1$hessian
  if (!is.null(H)) H <- 0.5 * (H + t(H))

  fittedparams <- c(
    R0 = exp(fit1$par[1]), v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]), c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    loglik = -fit1$value, convergence = fit1$convergence
  )

  list(parms = fittedparams, trans_parms = fit1$par, trans_hessian = H)
}

fit4_seir_1epic_poisson.loglikwithNPI_lognormal <- function(dat, init_het) {
  start_par <- c(log(2), log(1.1), log(12), logit(0.2))

  fit1 <- optim(
    par = start_par,
    fn  = f4_optim_1epic_poisson.loglikwithNPI_lognormal,
    sim.data = dat, init_het = init_het,
    method  = "Nelder-Mead",
    control = list(trace = 0, maxit = 1400),
    hessian = TRUE
  )

  H <- fit1$hessian
  if (!is.null(H)) H <- 0.5 * (H + t(H))

  fittedparams <- c(
    R0 = exp(fit1$par[1]), v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]), c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    loglik = -fit1$value, convergence = fit1$convergence
  )

  list(parms = fittedparams, trans_parms = fit1$par, trans_hessian = H)
}


# =============================================================================
# 9.  TWO-EPIDEMIC OBJECTIVE FUNCTIONS
# =============================================================================

f4_optim_2epic_poisson.loglikwithNPI <- function(par, sim.data_1, sim.data_2,
                                                  init_het_1, init_het_2) {
  z_R0 <- exp(par[1]);  z_v <- exp(par[2])
  z_t0 <- exp(par[3]);  z_c2 <- expit(par[4])

  alpha_fit <- 1 / (z_v^2)
  xval <- Discretize_gamma_LA_SM(
    n_groups = K, alpha = alpha_fit, beta = alpha_fit,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  if (xval$solution_text != "log_affine") return(.Machine$double.xmax)

  params <- c(xval$x,
              K = K, R0 = z_R0, gamma = gamma_spec, rho = rho_spec,
              delta = delta_spec, N = N,
              c_value1 = c_value1_spec, c_value2 = z_c2, c_value3 = c_value3_spec,
              t0 = z_t0, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec)
  names(params) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  ll1 <- poisson.loglik(params, sim.data = sim.data_1, init_het = init_het_1, times = sim.data_1$time)
  ll2 <- poisson.loglik(params, sim.data = sim.data_2, init_het = init_het_2, times = sim.data_2$time)
  -(ll1 + ll2)
}

f4_optim_2epic_poisson.loglikwithNPI_lognormal <- function(par, sim.data_1, sim.data_2,
                                                            init_het_1, init_het_2) {
  z_R0 <- exp(par[1]);  z_cv <- exp(par[2])
  z_t0 <- exp(par[3]);  z_c2 <- expit(par[4])
  if (!all(is.finite(c(z_R0, z_cv, z_t0, z_c2)))) return(.Machine$double.xmax)

  xval <- Discretize_lognormal_LA_SM(
    n_groups = K, CV = z_cv,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  if (is.null(xval$x) || length(xval$x) != K ||
      any(!is.finite(xval$x)) || any(xval$x <= 0)) return(.Machine$double.xmax)

  params <- c(xval$x,
              K = K, R0 = z_R0, gamma = gamma_spec, rho = rho_spec,
              delta = delta_spec, N = N,
              c_value1 = c_value1_spec, c_value2 = z_c2, c_value3 = c_value3_spec,
              t0 = z_t0, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec)
  names(params) <- c(paste0("x", 1:K), "K", "R0", "gamma", "rho", "delta", "N",
                     "c_value1", "c_value2", "c_value3", "t0", "t1", "t2", "t3")

  nll <- tryCatch({
    ll1 <- poisson.loglik.lognormal(params, sim.data = sim.data_1, init_het = init_het_1, times = sim.data_1$time)
    ll2 <- poisson.loglik.lognormal(params, sim.data = sim.data_2, init_het = init_het_2, times = sim.data_2$time)
    -(ll1 + ll2)
  }, error = function(e) .Machine$double.xmax)
  nll
}


# =============================================================================
# 10.  TWO-EPIDEMIC FITTING WRAPPERS
# =============================================================================

fit4_seir_2epic_poisson.loglik <- function(dat1, dat2, init_het_1, init_het_2) {
  start_par <- c(log(2), log(1.1), log(12), logit(0.2))

  fit1 <- optim(
    par = start_par,
    fn  = f4_optim_2epic_poisson.loglikwithNPI,
    sim.data_1 = dat1, sim.data_2 = dat2,
    init_het_1 = init_het_1, init_het_2 = init_het_2,
    method  = "Nelder-Mead",
    control = list(trace = 0, maxit = 1400),
    hessian = TRUE
  )

  H <- fit1$hessian
  if (!is.null(H)) H <- 0.5 * (H + t(H))

  fittedparams <- c(
    R0 = exp(fit1$par[1]), v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]), c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    loglik = -fit1$value, convergence = fit1$convergence
  )

  list(parms = fittedparams, trans_parms = fit1$par, trans_hessian = H)
}

fit4_seir_2epic_poisson.loglik_lognormal <- function(dat1, dat2, init_het_1, init_het_2) {
  start_par <- c(log(2), log(1.1), log(12), logit(0.2))

  fit1 <- optim(
    par = start_par,
    fn  = f4_optim_2epic_poisson.loglikwithNPI_lognormal,
    sim.data_1 = dat1, sim.data_2 = dat2,
    init_het_1 = init_het_1, init_het_2 = init_het_2,
    method  = "Nelder-Mead",
    control = list(trace = 0, maxit = 1400),
    hessian = TRUE
  )

  H <- fit1$hessian
  if (!is.null(H)) H <- 0.5 * (H + t(H))

  fittedparams <- c(
    R0 = exp(fit1$par[1]), v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]), c_value2 = expit(fit1$par[4]),
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    loglik = -fit1$value, convergence = fit1$convergence
  )

  list(parms = fittedparams, trans_parms = fit1$par, trans_hessian = H)
}

# =============================================================================
# 5.  HOMOGENEOUS REDUCED SEIR HELPERS
#
# These functions are the required homogeneous-model components from
# MLE_functions_paper.R used by the single- and two-epidemic misspecification
# scripts. They are kept here so the GitHub baseline analysis requires only one
# helper file.
# =============================================================================

Reduced.m_intervene <- function(t, y, parms) {
  
  with(as.list(c(t,y,parms)),{
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]
    
    # Piecewise-linear transmission multiplier c(t).
    prox <- c_value1
    
    if (t <= t0) {
      prox <- c_value1
    } else if (t <= t1) {
      prox <- c_value1 - (c_value1 - c_value2) * (t - t0) / (t1 - t0)
    } else if (t <= t2) {
      prox <- c_value2
    } else if (t<= t3) {
      prox <- c_value2 + (c_value3 - c_value2) * (t - t2) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    
    
    
    
    Beta <- R0*prox/(rho / delta + 1 / gamma)
    
    
    # SEIR model equations
    dS<- -Beta*(rho*E+I)*(S/N)^(1+v^2)    
    dE<- Beta*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI<- delta*E-gamma*I
    dR<- gamma*I
    dC<-delta*E
    return(list(c(dS, dE, dI, dR,dC)))
  })
  
}

# Time-aligned Poisson likelihood for the reduced SEIR model.
# Handles alignment between the model integration grid and the observation times.
poisson.loglik.withNPI.reduced <- function(params, sim.data, initial_state) {
  # Extract observation times
  times_data <- sim.data$time
  
  # Add time zero for ODE integration when observations start after zero
  if (min(times_data) > 0) {
    # Data starts at time 1, so add time 0 for ODE integration
    times_integration <- c(0, times_data)
    needs_trimming <- TRUE
  } else {
    # Data already includes time 0
    times_integration <- times_data
    needs_trimming <- FALSE
  }
  
  # Integrate the model equations
  out <- as.data.frame(ode(
    y = initial_state,
    times = times_integration,
    func = Reduced.m_intervene,
    parms = params
  ))
  
  # Remove the integration-only time-zero row if it was added
  if (needs_trimming) {
    out <- out[-1, ]  # Remove first row (time 0)
  }
  
  # Calculate daily incidence from cumulative cases
  # Use cumulative incidence differences aligned with the observation grid
  Daily_incidence <- diff(c(0, out[, "C"]))
  
  # Add incidence to dataframe
  df <- out %>% mutate(Inc = Daily_incidence)
  
  # Ensure values are valid for Poisson likelihood
  lambda_ <- ifelse(df[,"Inc"] <= 0, 0.0001, df[,"Inc"])
  
  # Validate sim.data
  if (!is.data.frame(sim.data)) {
    sim.data <- as.data.frame(sim.data)
  }
  
  # Confirm the required observation column is present
  if (!"reports" %in% names(sim.data)) {
    stop("The 'reports' column is missing in sim.data")
  }
  
  # Confirm that observations and model predictions are aligned
  if (nrow(sim.data) != length(lambda_)) {
    cat("Diagnostic: sim.data rows =", nrow(sim.data), ", lambda length =", length(lambda_), "\n")
    cat("Diagnostic: df rows =", nrow(df), "\n")
    cat("Diagnostic: sim.data time range =", range(sim.data$time), "\n")
    cat("Diagnostic: df time range =", range(df$time), "\n")
    stop("The lengths of sim.data and df do not match")
  }
  
  # Compute log likelihood
  loglik <- sum(dpois(
    x = sim.data[,"reports"],
    lambda = lambda_,
    log = TRUE
  ))
  
  return(loglik)
}

# ===============================================================================
# Homogeneous model functions (single epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f_optim_reducedm.poisloglikwithNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate the negative log-likelihood
  loglik <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  return(loglik)
}

#' Function to fit homogeneous model to a single epidemic with NPI
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_1epic_loglikwithNPI <- function(dat) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.3)), 
    fn = f_optim_reducedm.poisloglikwithNPI,
    sim.data = dat,
    method = "Nelder-Mead", 
    control = list(trace = 0, maxit = 1500),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

# ===============================================================================
# Homogeneous model functions (dual epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI for two epidemics
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data_1 First epidemic dataset
#' @param sim.data_2 Second epidemic dataset
#' @return Combined negative log-likelihood
f3_optim_reducedm.poisloglikwithNPI <- function(par, sim.data_1, sim.data_2) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate log-likelihood for first epidemic
  loglik1 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_1,
    initial_state = initial_state_1
  )
  
  # Calculate log-likelihood for second epidemic
  loglik2 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_2,
    initial_state = initial_state_2
  )
  
  # Sum negative log-likelihoods
  loglik = loglik1 + loglik2
  
  return(loglik)
}

#' Function to fit homogeneous model to two epidemics with NPI
#'
#' @param dat1 First epidemic dataset
#' @param dat2 Second epidemic dataset
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_2epic_loglikwithNPI <- function(dat1, dat2) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.2)), 
    fn = f3_optim_reducedm.poisloglikwithNPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 1600),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

