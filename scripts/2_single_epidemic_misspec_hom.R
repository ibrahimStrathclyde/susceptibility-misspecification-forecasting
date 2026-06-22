# =============================================================================
# SINGLE-EPIDEMIC MISSPECIFICATION ANALYSIS  (LA_SM version)
#
# =============================================================================

library(deSolve)
library(tidyverse)
library(MASS)
library(gridExtra)
library(grid)

# --- Reproducibility ---
SEED <- 20251217
set.seed(SEED)

# Source the updated helper functions (LA_SM version)
source("R/1_helpers.R")

# Ensure helpers exist
if (!exists("logit"))  logit  <- function(p) log(p / (1 - p))
if (!exists("expit"))  expit  <- function(x) 1 / (1 + exp(-x))

cat("=== SINGLE-EPIDEMIC MISSPECIFICATION ANALYSIS (LA_SM version) ===\n")

# =============================================================================
# STEP 1: FIXED PARAMETERS AND SETUP
# =============================================================================

N <- 100000
K <- 60
n_replicates <- 50

# Epidemiological parameters
R0_spec    <- 3.0
delta_spec <- 1 / 5.5
rho_spec   <- 0.5
gamma_spec <- 1 / 4

# Intervention parameters
t0_spec        <- 15
t1_spec        <- 20
t2_spec        <- 99
t3_spec        <- 100
tfinal_spec    <- 100L
c_value1_spec  <- 1
c_value2_spec  <- 0.3
c_value3_spec  <- 1

times_full <- 0:tfinal_spec

# Initial conditions
E0 <- 100
I0 <- 40

# Choose simulation distribution and CV value
SIMULATE_WITH <- "lognormal"      # Options: "gamma" or "lognormal"
CV_value      <- 1.414

cat("\nSimulation parameters:\n")
cat("- Simulating with:", SIMULATE_WITH, "distribution\n")
cat("- CV:", round(CV_value, 3), "\n")
cat("- R0:", R0_spec, "\n")
cat("- NPI strength:", c_value2_spec, "\n")
cat("- Number of replicates:", n_replicates, "\n")

# =============================================================================
# STEP 1b: GLOBALS REQUIRED BY THE HOMOGENEOUS FIT
#
#  fit3_hom_1epic_loglikwithNPI() calls f_optim_reducedm.poisloglikwithNPI,
#  which looks up v_spec, initial_state, and the t*_spec / c_value*_spec /
#  rate parameters in .GlobalEnv. We set v_spec = 0 for the homogeneous case
#  and build initial_state from the single-epidemic initial conditions.
#
#  This block must run BEFORE Step 3, otherwise the homogeneous fit will fail
#  with "object 'v_spec' not found".
# =============================================================================

v_spec <- 0
assign("v_spec", v_spec, envir = .GlobalEnv)

initial_state <- c(S = N - E0 - I0, E = E0, I = I0, R = 0, C = 0)
assign("initial_state", initial_state, envir = .GlobalEnv)

# Re-export everything the homogeneous objective function reads from .GlobalEnv.
# (These are visible here, but being explicit protects against sourcing order
#  surprises when the script is wrapped in a function or a createfile.)
for (const in c("N", "t0_spec", "t1_spec", "t2_spec", "t3_spec",
              "c_value1_spec", "c_value2_spec", "c_value3_spec",
              "rho_spec", "delta_spec", "gamma_spec", "tfinal_spec")) {
  assign(const, get(const), envir = .GlobalEnv)
}
rm(const)

# =============================================================================
# STEP 2: SIMULATE DATA
# =============================================================================

cat("\n=== SIMULATING DATA ===\n")

simulated_datasets <- list()

pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
for (i in 1:n_replicates) {
  setTxtProgressBar(pb, i)
  
  if (SIMULATE_WITH == "gamma") {
    sim_result <- simulate_cases_hetsus_model_gamma(
      alpha = 1 / (CV_value^2), K = K,
      R0 = R0_spec, delta = delta_spec, rho = rho_spec, gamma = gamma_spec,
      N = N, E0 = E0, I0 = I0,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
  } else {
    sim_result <- simulate_cases_hetsus_model_lognormal(
      CV = CV_value, K = K,
      R0 = R0_spec, delta = delta_spec, rho = rho_spec, gamma = gamma_spec,
      N = N, E0 = E0, I0 = I0,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
  }
  
  sim_data <- sim_result$sim_data
  names(sim_data) <- c("time", "reports")
  simulated_datasets[[i]] <- sim_data
}
close(pb)
cat("\nData simulation complete!\n")

# =============================================================================
# STEP 3: FIT BOTH GAMMA AND LOGNORMAL MODELS
# =============================================================================

cat("\n=== FITTING MODELS ===\n")

results_gamma    <- NULL
results_lognormal <- NULL
results_homogeneous <- NULL
gamma_covariances    <- list()
lognormal_covariances <- list()
hom_covariances       <- list()

pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
for (i in 1:n_replicates) {
  setTxtProgressBar(pb, i)
  current_data <- simulated_datasets[[i]]
  
  # ---------- FIT GAMMA ----------
  xval_gam <- Discretize_gamma_LA_SM(
    n_groups = K,
    alpha = 1 / (CV_value^2), beta = 1 / (CV_value^2),
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  init_het_gamma <- initDistr(q = xval_gam$q, K = K, E0 = E0, I0 = I0, N = N)
  
  z_mle_gamma <- tryCatch({
    fit4_seir_1epic_poisson.loglikwithNPI(dat = current_data, init_het = init_het_gamma)
  }, error = function(e) NULL)
  
  if (!is.null(z_mle_gamma) && z_mle_gamma$parms["convergence"] == 0) {
    gamma_fit <- data.frame(
      dataset_id = i,
      R0 = z_mle_gamma$parms["R0"], v = z_mle_gamma$parms["v"],
      t0 = z_mle_gamma$parms["t0"], c_value2 = z_mle_gamma$parms["c_value2"],
      aic = z_mle_gamma$parms["AIC"], convergence = z_mle_gamma$parms["convergence"],
      fit_dist = "gamma", stringsAsFactors = FALSE
    )
    gamma_fit$hess_pd <- FALSE;  gamma_fit$condition_number <- NA
    
    if (!is.null(z_mle_gamma$trans_hessian)) {
      tryCatch({
        H  <- z_mle_gamma$trans_hessian
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        tol <- .Machine$double.eps^(2 / 3)
        if (all(ev > tol)) {
          gamma_fit$hess_pd <- TRUE
          gamma_fit$condition_number <- max(ev) / min(ev)
          V <- chol2inv(chol(H));  z_se <- sqrt(diag(V))
          z_corr_matrix <- cov2cor(V)
          gamma_covariances[[i]] <- V
          
          par_ucl <- z_mle_gamma$trans_parms + 1.96 * z_se
          par_lcl <- z_mle_gamma$trans_parms - 1.96 * z_se
          
          gamma_fit$R0_lcl       <- exp(par_lcl[1]);   gamma_fit$R0_ucl       <- exp(par_ucl[1])
          gamma_fit$v_lcl        <- exp(par_lcl[2]);   gamma_fit$v_ucl        <- exp(par_ucl[2])
          gamma_fit$t0_lcl       <- exp(par_lcl[3]);   gamma_fit$t0_ucl       <- exp(par_ucl[3])
          gamma_fit$c_value2_lcl <- expit(par_lcl[4]); gamma_fit$c_value2_ucl <- expit(par_ucl[4])
          gamma_fit$v_c_corr     <- z_corr_matrix[2, 4]
        }
      }, error = function(e) NULL)
    }
    results_gamma <- if (is.null(results_gamma)) gamma_fit else bind_rows(results_gamma, gamma_fit)
  }
  
  # ---------- FIT LOGNORMAL ----------
  val <- Discretize_lognormal_LA_SM(
    n_groups = K, CV = CV_value,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  init_het_lognorm <- initDistr(q = val$q, K = K, E0 = E0, I0 = I0, N = N)
  
  z_mle_lognormal <- tryCatch({
    fit4_seir_1epic_poisson.loglikwithNPI_lognormal(dat = current_data, init_het = init_het_lognorm)
  }, error = function(e) NULL)
  
  if (!is.null(z_mle_lognormal) && z_mle_lognormal$parms["convergence"] == 0) {
    lognormal_fit <- data.frame(
      dataset_id = i,
      R0 = z_mle_lognormal$parms["R0"], v = z_mle_lognormal$parms["v"],
      t0 = z_mle_lognormal$parms["t0"], c_value2 = z_mle_lognormal$parms["c_value2"],
      aic = z_mle_lognormal$parms["AIC"], convergence = z_mle_lognormal$parms["convergence"],
      fit_dist = "lognormal", stringsAsFactors = FALSE
    )
    lognormal_fit$hess_pd <- FALSE;  lognormal_fit$condition_number <- NA
    
    if (!is.null(z_mle_lognormal$trans_hessian)) {
      tryCatch({
        H  <- z_mle_lognormal$trans_hessian
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        tol <- .Machine$double.eps^(2 / 3)
        if (all(ev > tol)) {
          lognormal_fit$hess_pd <- TRUE
          lognormal_fit$condition_number <- max(ev) / min(ev)
          V <- chol2inv(chol(H));  z_se <- sqrt(diag(V))
          z_corr_matrix <- cov2cor(V)
          lognormal_covariances[[i]] <- V
          
          par_ucl <- z_mle_lognormal$trans_parms + 1.96 * z_se
          par_lcl <- z_mle_lognormal$trans_parms - 1.96 * z_se
          
          lognormal_fit$R0_lcl       <- exp(par_lcl[1]);   lognormal_fit$R0_ucl       <- exp(par_ucl[1])
          lognormal_fit$v_lcl        <- exp(par_lcl[2]);   lognormal_fit$v_ucl        <- exp(par_ucl[2])
          lognormal_fit$t0_lcl       <- exp(par_lcl[3]);   lognormal_fit$t0_ucl       <- exp(par_ucl[3])
          lognormal_fit$c_value2_lcl <- expit(par_lcl[4]); lognormal_fit$c_value2_ucl <- expit(par_ucl[4])
          lognormal_fit$v_c_corr     <- z_corr_matrix[2, 4]
        }
      }, error = function(e) NULL)
    }
    results_lognormal <- if (is.null(results_lognormal)) lognormal_fit else bind_rows(results_lognormal, lognormal_fit)
  }
  
  # ---------- FIT HOMOGENEOUS (reference baseline) ----------
  # No heterogeneity: v is fixed at 0 (set in .GlobalEnv above) and the
  # objective function uses Reduced.m_intervene over the 5-state SEIR system.
  # This fit has 3 estimated parameters (R0, t0, c_value2) instead of 4.
  z_mle_hom <- tryCatch({
    fit3_hom_1epic_loglikwithNPI(dat = current_data)
  }, error = function(e) NULL)
  
  if (!is.null(z_mle_hom) && z_mle_hom$parms["convergence"] == 0) {
    hom_fit <- data.frame(
      dataset_id = i,
      R0 = z_mle_hom$parms["R0"], v = NA_real_,
      t0 = z_mle_hom$parms["t0"], c_value2 = z_mle_hom$parms["c_value2"],
      aic = z_mle_hom$parms["AIC"], convergence = z_mle_hom$parms["convergence"],
      fit_dist = "homogeneous", stringsAsFactors = FALSE
    )
    hom_fit$hess_pd <- FALSE;  hom_fit$condition_number <- NA
    
    if (!is.null(z_mle_hom$trans_hessian)) {
      tryCatch({
        H  <- z_mle_hom$trans_hessian
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        tol <- .Machine$double.eps^(2 / 3)
        if (all(ev > tol)) {
          hom_fit$hess_pd <- TRUE
          hom_fit$condition_number <- max(ev) / min(ev)
          V <- chol2inv(chol(H));  z_se <- sqrt(diag(V))
          hom_covariances[[i]] <- V
          
          par_ucl <- z_mle_hom$trans_parms + 1.96 * z_se
          par_lcl <- z_mle_hom$trans_parms - 1.96 * z_se
          
          # The homogeneous fit has only 3 estimated parameters, so the
          # trans_parms vector is (log R0, log t0, logit c_value2).
          # v has no interval (fixed at 0) and there is no v-c correlation.
          hom_fit$R0_lcl       <- exp(par_lcl[1]);   hom_fit$R0_ucl       <- exp(par_ucl[1])
          hom_fit$v_lcl        <- NA_real_;          hom_fit$v_ucl        <- NA_real_
          hom_fit$t0_lcl       <- exp(par_lcl[2]);   hom_fit$t0_ucl       <- exp(par_ucl[2])
          hom_fit$c_value2_lcl <- expit(par_lcl[3]); hom_fit$c_value2_ucl <- expit(par_ucl[3])
          hom_fit$v_c_corr     <- NA_real_
        }
      }, error = function(e) NULL)
    }
    results_homogeneous <- if (is.null(results_homogeneous)) hom_fit else bind_rows(results_homogeneous, hom_fit)
  }
}
close(pb)
cat("\nModel fitting complete!\n")

# =============================================================================
# STEP 4: ANALYZE RESULTS
# =============================================================================

cat("\n=== ANALYZING RESULTS ===\n")

valid_gamma       <- results_gamma       %>% filter(hess_pd == TRUE, convergence == 0)
valid_lognormal   <- results_lognormal   %>% filter(hess_pd == TRUE, convergence == 0)
valid_homogeneous <- if (!is.null(results_homogeneous)) {
  results_homogeneous %>% filter(hess_pd == TRUE, convergence == 0)
} else {
  results_homogeneous[0, , drop = FALSE]  # empty frame with the right schema
}

cat("\nConvergence summary:\n")
cat("- Gamma model: ",       nrow(valid_gamma),       "out of", n_replicates, "converged (PD Hessian)\n")
cat("- Lognormal model: ",   nrow(valid_lognormal),   "out of", n_replicates, "converged (PD Hessian)\n")
cat("- Homogeneous model: ", nrow(valid_homogeneous), "out of", n_replicates, "converged (PD Hessian)\n")





# =============================================================================
# BLOCK 1: Extract parameter summary table from valid_* objects in memory
# Works for both single-epidemic and two-epidemic scripts
# =============================================================================

# True values
# True values
R0_true <- 3.0;  CV_true <- 1.414;  t0_true <- 15.0;  c2_true <- 0.3

summarise_model <- function(df, label) {
  
  safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
  safe_sd     <- function(x) if (all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE)
  safe_mean   <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  safe_cov    <- function(lcl, ucl, truth) {
    ok <- !is.na(lcl) & !is.na(ucl)
    if (sum(ok) == 0) return(NA_real_)
    mean(lcl[ok] <= truth & truth <= ucl[ok]) * 100
  }
  
  data.frame(
    Model     = label,
    Parameter = c("R0", "nu", "t0", "c_star"),
    
    Median = c(safe_median(df$R0), safe_median(df$v),
               safe_median(df$t0), safe_median(df$c_value2)),
    
    SD     = c(safe_sd(df$R0), safe_sd(df$v),
               safe_sd(df$t0), safe_sd(df$c_value2)),
    
    CI_Width = c(safe_mean(df$R0_ucl - df$R0_lcl),
                 safe_mean(df$v_ucl  - df$v_lcl),
                 safe_mean(df$t0_ucl - df$t0_lcl),
                 safe_mean(df$c_value2_ucl - df$c_value2_lcl)),
    
    Coverage = c(safe_cov(df$R0_lcl, df$R0_ucl, R0_true),
                 safe_cov(df$v_lcl,  df$v_ucl,  CV_true),
                 safe_cov(df$t0_lcl, df$t0_ucl, t0_true),
                 safe_cov(df$c_value2_lcl, df$c_value2_ucl, c2_true)),
    
    stringsAsFactors = FALSE
  )
}

# --- Single epidemic ---
tab_gamma_1epi     <- summarise_model(valid_gamma,       "Gamma")
tab_lognormal_1epi <- summarise_model(valid_lognormal,   "Lognormal")
tab_hom_1epi       <- summarise_model(valid_homogeneous, "Homogeneous")

summary_1epi <- bind_rows(tab_gamma_1epi, tab_lognormal_1epi, tab_hom_1epi)

# Format with NA-safe sprintf
summary_1epi$Median_SD <- ifelse(
  is.na(summary_1epi$Median),
  "---",
  sprintf("%.3f (%.3f)", summary_1epi$Median, summary_1epi$SD)
)
summary_1epi$Width_fmt <- ifelse(
  is.na(summary_1epi$CI_Width),
  "---",
  sprintf("%.3f", summary_1epi$CI_Width)
)
summary_1epi$Cov_fmt <- ifelse(
  is.na(summary_1epi$Coverage),
  "---",
  sprintf("%.1f", summary_1epi$Coverage)
)

cat("\n=== SINGLE EPIDEMIC PARAMETER SUMMARY ===\n")
print(summary_1epi %>% dplyr::select(Model, Parameter, Median_SD, Width_fmt, Cov_fmt))



summary_gamma <- valid_gamma %>%
  summarise(
    Model = "Gamma",
    R0_mean = median(R0, na.rm = TRUE),
    R0_bias = median(R0 - R0_spec, na.rm = TRUE),
    R0_rel_bias = median((R0 - R0_spec) / R0_spec * 100, na.rm = TRUE),
    v_mean = median(v, na.rm = TRUE),
    v_bias = median(v - CV_value, na.rm = TRUE),
    v_rel_bias = median((v - CV_value) / CV_value * 100, na.rm = TRUE),
    t0_mean = median(t0, na.rm = TRUE),
    t0_bias = median(t0 - t0_spec, na.rm = TRUE),
    c_value2_mean = median(c_value2, na.rm = TRUE),
    c_value2_bias = median(c_value2 - c_value2_spec, na.rm = TRUE),
    c_value2_rel_bias = median((c_value2 - c_value2_spec) / c_value2_spec * 100, na.rm = TRUE),
    R0_coverage       = mean(R0_lcl <= R0_spec & R0_spec <= R0_ucl, na.rm = TRUE),
    v_coverage        = mean(v_lcl  <= CV_value & CV_value <= v_ucl, na.rm = TRUE),
    t0_coverage       = mean(t0_lcl <= t0_spec  & t0_spec  <= t0_ucl, na.rm = TRUE),
    c_value2_coverage = mean(c_value2_lcl <= c_value2_spec & c_value2_spec <= c_value2_ucl, na.rm = TRUE),
    median_condition  = median(condition_number, na.rm = TRUE),
    median_v_c_corr   = median(v_c_corr, na.rm = TRUE),
    n_converged       = n()
  )

summary_lognormal <- valid_lognormal %>%
  summarise(
    Model = "Lognormal",
    R0_mean = median(R0, na.rm = TRUE),
    R0_bias = median(R0 - R0_spec, na.rm = TRUE),
    R0_rel_bias = median((R0 - R0_spec) / R0_spec * 100, na.rm = TRUE),
    v_mean = median(v, na.rm = TRUE),
    v_bias = median(v - CV_value, na.rm = TRUE),
    v_rel_bias = median((v - CV_value) / CV_value * 100, na.rm = TRUE),
    t0_mean = median(t0, na.rm = TRUE),
    t0_bias = median(t0 - t0_spec, na.rm = TRUE),
    c_value2_mean = median(c_value2, na.rm = TRUE),
    c_value2_bias = median(c_value2 - c_value2_spec, na.rm = TRUE),
    c_value2_rel_bias = median((c_value2 - c_value2_spec) / c_value2_spec * 100, na.rm = TRUE),
    R0_coverage       = mean(R0_lcl <= R0_spec & R0_spec <= R0_ucl, na.rm = TRUE),
    v_coverage        = mean(v_lcl  <= CV_value & CV_value <= v_ucl, na.rm = TRUE),
    t0_coverage       = mean(t0_lcl <= t0_spec  & t0_spec  <= t0_ucl, na.rm = TRUE),
    c_value2_coverage = mean(c_value2_lcl <= c_value2_spec & c_value2_spec <= c_value2_ucl, na.rm = TRUE),
    median_condition  = median(condition_number, na.rm = TRUE),
    median_v_c_corr   = median(v_c_corr, na.rm = TRUE),
    n_converged       = n()
  )

summary_table <- bind_rows(summary_gamma, summary_lognormal)



#~RELATIVE BIAS (%):
  print(summary_table %>% dplyr::select(Model, R0_rel_bias, v_rel_bias, c_value2_rel_bias))

# --- Export combined valid results for external plotting ---
valid_gamma$Model       <- "Gamma"
valid_lognormal$Model   <- "Lognormal"
valid_homogeneous$Model <- "Homogeneous"

combined_data <- bind_rows(valid_gamma, valid_lognormal, valid_homogeneous)
combined_data$Model <- factor(combined_data$Model, levels = c("Gamma", "Lognormal", "Homogeneous"))

write.csv(combined_data,
          paste0("combined_valid_single_epi_sim_", SIMULATE_WITH, ".csv"),
          row.names = FALSE)
cat("Saved: combined_valid_single_epi_sim_", SIMULATE_WITH, ".csv\n", sep = "")


# --- Homogeneous summary (v-related columns are NA because v is fixed at 0) ---
if (nrow(valid_homogeneous) > 0) {
  summary_homogeneous <- valid_homogeneous %>%
    summarise(
      Model = "Homogeneous",
      R0_mean = median(R0, na.rm = TRUE),
      R0_bias = median(R0 - R0_spec, na.rm = TRUE),
      R0_rel_bias = median((R0 - R0_spec) / R0_spec * 100, na.rm = TRUE),
      v_mean = NA_real_, v_bias = NA_real_, v_rel_bias = NA_real_,
      t0_mean = median(t0, na.rm = TRUE),
      t0_bias = median(t0 - t0_spec, na.rm = TRUE),
      c_value2_mean = median(c_value2, na.rm = TRUE),
      c_value2_bias = median(c_value2 - c_value2_spec, na.rm = TRUE),
      c_value2_rel_bias = median((c_value2 - c_value2_spec) / c_value2_spec * 100, na.rm = TRUE),
      R0_coverage       = mean(R0_lcl <= R0_spec & R0_spec <= R0_ucl, na.rm = TRUE),
      v_coverage        = NA_real_,
      t0_coverage       = mean(t0_lcl <= t0_spec  & t0_spec  <= t0_ucl, na.rm = TRUE),
      c_value2_coverage = mean(c_value2_lcl <= c_value2_spec & c_value2_spec <= c_value2_ucl, na.rm = TRUE),
      median_condition  = median(condition_number, na.rm = TRUE),
      median_v_c_corr   = NA_real_,
      n_converged       = n()
    )
  summary_table <- bind_rows(summary_table, summary_homogeneous)
}

cat("\nSUMMARY RESULTS:\n")
cat("Data simulated with:", SIMULATE_WITH, "distribution (CV =", CV_value, ")\n\n")
print(summary_table %>% dplyr::select(Model, R0_mean, v_mean, t0_mean, c_value2_mean))
cat("\nRELATIVE BIAS (%):\n")
print(summary_table %>% dplyr::select(Model, R0_rel_bias, v_rel_bias, c_value2_rel_bias))
cat("\nCOVERAGE:\n")
print(summary_table %>% dplyr::select(Model, R0_coverage, v_coverage, t0_coverage, c_value2_coverage))
cat("\nIDENTIFIABILITY:\n")
print(summary_table %>% dplyr::select(Model, median_condition, median_v_c_corr))

# Save CSV results
write.csv(results_gamma,     paste0("results_gamma_sim_",     SIMULATE_WITH, ".csv"), row.names = FALSE)
write.csv(results_lognormal, paste0("results_lognormal_sim_", SIMULATE_WITH, ".csv"), row.names = FALSE)
if (!is.null(results_homogeneous)) {
  write.csv(results_homogeneous, paste0("results_homogeneous_sim_", SIMULATE_WITH, ".csv"), row.names = FALSE)
}
write.csv(summary_table,     paste0("summary_sim_",           SIMULATE_WITH, ".csv"), row.names = FALSE)
saveRDS(gamma_covariances,    paste0("gamma_covariances_sim_",    SIMULATE_WITH, ".rds"))
saveRDS(lognormal_covariances, paste0("lognormal_covariances_sim_", SIMULATE_WITH, ".rds"))
saveRDS(hom_covariances,       paste0("hom_covariances_sim_",       SIMULATE_WITH, ".rds"))
# =============================================================================
# STEP 5: PREDICTION TRAJECTORIES AND FORECAST METRICS
# =============================================================================

T_fit            <- 100L
tfinal_forecast  <- 250L
times_full       <- 0:tfinal_forecast
times_eval       <- (T_fit + 1):tfinal_forecast
N_SAMPLES        <- 300L

TRUTH_DIST <- toupper(SIMULATE_WITH)
TRUTH_CV   <- CV_value
sw          <- SIMULATE_WITH

# -----------------------------------------------------------------------------
# 5.1 Observed training data: mean across simulated replicates
# -----------------------------------------------------------------------------

obs_col <- {
  nm <- names(simulated_datasets[[1]])
  if ("reports" %in% nm) {
    "reports"
  } else if ("cases" %in% nm) {
    "cases"
  } else {
    NA_character_
  }
}

if (is.na(obs_col)) {
  stop("The simulated datasets must contain either a 'reports' or 'cases' column.")
}

dat_fit <- dplyr::bind_rows(simulated_datasets) %>%
  dplyr::filter(time >= 1, time <= T_fit) %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(
    reports = mean(.data[[obs_col]], na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 5.2 ODE output helpers
# -----------------------------------------------------------------------------

daily_incidence_het <- function(params_vec, init_state, times, K) {
  out <- as.data.frame(deSolve::ode(
    y = init_state,
    times = times,
    func = hetsus_model.ct,
    parms = params_vec
  ))
  
  C_cols <- grep("^C[0-9]+$", names(out))
  if (length(C_cols) != K) {
    stop("Could not identify all cumulative-incidence columns C1, ..., CK.")
  }
  
  Ccum <- rowSums(out[, C_cols, drop = FALSE], na.rm = TRUE)
  
  tibble::tibble(
    time = out$time,
    Inc  = pmax(c(0, diff(Ccum)), 0)
  )
}

daily_incidence_hom <- function(params_vec, init_state, times) {
  out <- as.data.frame(deSolve::ode(
    y = init_state,
    times = times,
    func = Reduced.m_intervene,
    parms = params_vec
  ))
  
  tibble::tibble(
    time = out$time,
    Inc  = pmax(c(0, diff(out$C)), 0)
  )
}

# -----------------------------------------------------------------------------
# 5.3 Parameter-vector builders
# -----------------------------------------------------------------------------

create_gamma_components <- function(R0, v, t0, c2) {
  alpha <- 1 / (v^2)
  
  d <- Discretize_gamma_LA_SM(
    n_groups = K,
    alpha = alpha,
    beta = alpha,
    spacing = "equal",
    rep = "condexp",
    calibration = "log_affine",
    print_message = FALSE
  )
  
  params <- c(
    d$x,
    K = K,
    R0 = R0,
    gamma = gamma_spec,
    rho = rho_spec,
    delta = delta_spec,
    N = N,
    c_value1 = c_value1_spec,
    c_value2 = c2,
    c_value3 = c_value3_spec,
    t0 = t0,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec
  )
  
  names(params) <- c(
    paste0("x", seq_len(K)),
    "K", "R0", "gamma", "rho", "delta", "N",
    "c_value1", "c_value2", "c_value3",
    "t0", "t1", "t2", "t3"
  )
  
  list(params = params, q = d$q)
}

create_lognormal_components <- function(R0, v, t0, c2) {
  d <- Discretize_lognormal_LA_SM(
    n_groups = K,
    CV = v,
    spacing = "equal",
    rep = "condexp",
    calibration = "log_affine",
    print_message = FALSE
  )
  
  params <- c(
    d$x,
    K = K,
    R0 = R0,
    gamma = gamma_spec,
    rho = rho_spec,
    delta = delta_spec,
    N = N,
    c_value1 = c_value1_spec,
    c_value2 = c2,
    c_value3 = c_value3_spec,
    t0 = t0,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec
  )
  
  names(params) <- c(
    paste0("x", seq_len(K)),
    "K", "R0", "gamma", "rho", "delta", "N",
    "c_value1", "c_value2", "c_value3",
    "t0", "t1", "t2", "t3"
  )
  
  list(params = params, q = d$q)
}

create_hom_params <- function(R0, t0, c2) {
  c(
    R0 = R0,
    v = 0,
    t0 = t0,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c2,
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N,
    tfinal = tfinal_forecast
  )
}

# -----------------------------------------------------------------------------
# 5.4 Forecast-window summary metrics
# -----------------------------------------------------------------------------

compute_peak_metrics <- function(traj, window) {
  sub <- traj[traj$time %in% window, , drop = FALSE]
  
  if (nrow(sub) == 0 || all(!is.finite(sub$Inc))) {
    return(list(
      peak_height = NA_real_,
      peak_day    = NA_real_,
      final_size  = NA_real_
    ))
  }
  
  idx <- which.max(sub$Inc)
  
  list(
    peak_height = sub$Inc[idx],
    peak_day    = sub$time[idx],
    final_size  = sum(sub$Inc, na.rm = TRUE)
  )
}

percent_difference <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || b == 0) {
    return(NA_real_)
  }
  100 * (a - b) / b
}

# -----------------------------------------------------------------------------
# 5.5 Deterministic truth trajectory
# -----------------------------------------------------------------------------

truth_components <- if (SIMULATE_WITH == "gamma") {
  create_gamma_components(
    R0 = R0_spec,
    v = CV_value,
    t0 = t0_spec,
    c2 = c_value2_spec
  )
} else {
  create_lognormal_components(
    R0 = R0_spec,
    v = CV_value,
    t0 = t0_spec,
    c2 = c_value2_spec
  )
}

init_truth <- initDistr(
  q = truth_components$q,
  K = K,
  E0 = E0,
  I0 = I0,
  N = N
)

traj_true <- daily_incidence_het(
  params_vec = truth_components$params,
  init_state = init_truth,
  times = times_full,
  K = K
) %>%
  dplyr::filter(time >= 1)

truth_metrics <- compute_peak_metrics(traj_true, times_eval)

cat(sprintf(
  "\nTruth forecast peak: height = %.1f, day = %d, final_size = %.0f\n",
  truth_metrics$peak_height,
  truth_metrics$peak_day,
  truth_metrics$final_size
))

# -----------------------------------------------------------------------------
# 5.6 Per-replicate forecast metrics
# -----------------------------------------------------------------------------

forecast_peak_from_fit <- function(fit, fitted_model) {
  
  if (fitted_model == "gamma") {
    
    comp <- create_gamma_components(
      R0 = fit$R0,
      v = fit$v,
      t0 = fit$t0,
      c2 = fit$c_value2
    )
    
    init_fit <- initDistr(
      q = comp$q,
      K = K,
      E0 = E0,
      I0 = I0,
      N = N
    )
    
    traj_fit <- daily_incidence_het(
      params_vec = comp$params,
      init_state = init_fit,
      times = times_full,
      K = K
    )
    
  } else if (fitted_model == "lognormal") {
    
    comp <- create_lognormal_components(
      R0 = fit$R0,
      v = fit$v,
      t0 = fit$t0,
      c2 = fit$c_value2
    )
    
    init_fit <- initDistr(
      q = comp$q,
      K = K,
      E0 = E0,
      I0 = I0,
      N = N
    )
    
    traj_fit <- daily_incidence_het(
      params_vec = comp$params,
      init_state = init_fit,
      times = times_full,
      K = K
    )
    
  } else if (fitted_model == "homogeneous") {
    
    params_fit <- create_hom_params(
      R0 = fit$R0,
      t0 = fit$t0,
      c2 = fit$c_value2
    )
    
    traj_fit <- daily_incidence_hom(
      params_vec = params_fit,
      init_state = initial_state,
      times = times_full
    )
    
  } else {
    stop("Unknown fitted_model: ", fitted_model)
  }
  
  traj_fit <- traj_fit %>%
    dplyr::filter(time >= 1)
  
  compute_peak_metrics(traj_fit, times_eval)
}

create_peak_diff_row <- function(fit, fitted_model, metrics_fit) {
  data.frame(
    CV_true = CV_value,
    truth_dist = SIMULATE_WITH,
    fitted_model = fitted_model,
    dataset_id = fit$dataset_id,
    
    R0_fit = fit$R0,
    v_fit = fit$v,
    t0_fit = fit$t0,
    c2_fit = fit$c_value2,
    
    peak_height_fit = metrics_fit$peak_height,
    peak_day_fit = metrics_fit$peak_day,
    final_size_fit = metrics_fit$final_size,
    
    peak_height_truth = truth_metrics$peak_height,
    peak_day_truth = truth_metrics$peak_day,
    final_size_truth = truth_metrics$final_size,
    
    peak_height_rel_diff_pct = percent_difference(
      metrics_fit$peak_height,
      truth_metrics$peak_height
    ),
    
    peak_day_diff = if (
      is.finite(metrics_fit$peak_day) &&
      is.finite(truth_metrics$peak_day)
    ) {
      metrics_fit$peak_day - truth_metrics$peak_day
    } else {
      NA_real_
    },
    
    final_size_rel_diff_pct = percent_difference(
      metrics_fit$final_size,
      truth_metrics$final_size
    ),
    
    stringsAsFactors = FALSE
  )
}

process_model_fits <- function(valid_fits, fitted_model) {
  if (is.null(valid_fits) || nrow(valid_fits) == 0) {
    return(NULL)
  }
  
  rows <- vector("list", nrow(valid_fits))
  
  for (i in seq_len(nrow(valid_fits))) {
    fit_i <- valid_fits[i, ]
    
    metrics_i <- tryCatch(
      forecast_peak_from_fit(fit_i, fitted_model),
      error = function(e) NULL
    )
    
    if (!is.null(metrics_i)) {
      rows[[i]] <- create_peak_diff_row(fit_i, fitted_model, metrics_i)
    }
  }
  
  dplyr::bind_rows(rows)
}

cat("\n=== COMPUTING PER-REPLICATE FORECAST PEAK METRICS ===\n")

peak_diff_single <- dplyr::bind_rows(
  process_model_fits(valid_gamma, "gamma"),
  process_model_fits(valid_lognormal, "lognormal"),
  process_model_fits(valid_homogeneous, "homogeneous")
)

write.csv(
  peak_diff_single,
  paste0("peak_diff_single_sim_", SIMULATE_WITH, ".csv"),
  row.names = FALSE
)

cat(sprintf(
  "Wrote peak_diff_single_sim_%s.csv  (rows: %d)\n",
  SIMULATE_WITH,
  nrow(peak_diff_single)
))

cat("\nForecast-window medians relative to deterministic truth:\n")
print(
  aggregate(
    cbind(
      peak_height_rel_diff_pct,
      peak_day_diff,
      final_size_rel_diff_pct
    ) ~ fitted_model,
    data = peak_diff_single,
    FUN = function(x) median(x, na.rm = TRUE)
  )
)

# -----------------------------------------------------------------------------
# 5.7 Median parameter estimates for mean trajectories
# -----------------------------------------------------------------------------

clip_probability <- function(x) {
  pmin(pmax(x, 1e-8), 1 - 1e-8)
}

median_het_fit <- function(valid_fits, model_name) {
  if (is.null(valid_fits) || nrow(valid_fits) == 0) {
    stop("No valid ", model_name, " fits are available.")
  }
  
  list(
    R0 = max(stats::median(valid_fits$R0, na.rm = TRUE), 1e-8),
    v  = max(stats::median(valid_fits$v, na.rm = TRUE), 1e-8),
    t0 = max(stats::median(valid_fits$t0, na.rm = TRUE), 1e-8),
    c2 = clip_probability(stats::median(valid_fits$c_value2, na.rm = TRUE))
  )
}

median_hom_fit <- function(valid_fits) {
  if (is.null(valid_fits) || nrow(valid_fits) == 0) {
    return(NULL)
  }
  
  list(
    R0 = max(stats::median(valid_fits$R0, na.rm = TRUE), 1e-8),
    t0 = max(stats::median(valid_fits$t0, na.rm = TRUE), 1e-8),
    c2 = clip_probability(stats::median(valid_fits$c_value2, na.rm = TRUE))
  )
}

fit_gam <- median_het_fit(valid_gamma, "Gamma")
fit_log <- median_het_fit(valid_lognormal, "Lognormal")
fit_hom <- median_hom_fit(valid_homogeneous)

# -----------------------------------------------------------------------------
# 5.8 Mean fitted trajectories
# -----------------------------------------------------------------------------

gam_components <- create_gamma_components(
  R0 = fit_gam$R0,
  v = fit_gam$v,
  t0 = fit_gam$t0,
  c2 = fit_gam$c2
)

log_components <- create_lognormal_components(
  R0 = fit_log$R0,
  v = fit_log$v,
  t0 = fit_log$t0,
  c2 = fit_log$c2
)

init_gam <- initDistr(
  q = gam_components$q,
  K = K,
  E0 = E0,
  I0 = I0,
  N = N
)

init_log <- initDistr(
  q = log_components$q,
  K = K,
  E0 = E0,
  I0 = I0,
  N = N
)

traj_gam <- daily_incidence_het(
  params_vec = gam_components$params,
  init_state = init_gam,
  times = times_full,
  K = K
) %>%
  dplyr::filter(time >= 1)

traj_log <- daily_incidence_het(
  params_vec = log_components$params,
  init_state = init_log,
  times = times_full,
  K = K
) %>%
  dplyr::filter(time >= 1)

if (!is.null(fit_hom)) {
  params_hom <- create_hom_params(
    R0 = fit_hom$R0,
    t0 = fit_hom$t0,
    c2 = fit_hom$c2
  )
  
  traj_hom <- daily_incidence_hom(
    params_vec = params_hom,
    init_state = initial_state,
    times = times_full
  ) %>%
    dplyr::filter(time >= 1)
} else {
  traj_hom <- NULL
}

# -----------------------------------------------------------------------------
# 5.9 Forecast uncertainty bands
# -----------------------------------------------------------------------------

median_covariance <- function(cov_list, dimension) {
  mats <- Filter(
    function(M) {
      is.matrix(M) &&
        all(dim(M) == c(dimension, dimension)) &&
        all(is.finite(M))
    },
    cov_list
  )
  
  if (!length(mats)) {
    return(NULL)
  }
  
  med <- matrix(NA_real_, dimension, dimension)
  
  for (i in seq_len(dimension)) {
    for (j in seq_len(dimension)) {
      med[i, j] <- stats::median(
        vapply(mats, function(M) M[i, j], numeric(1)),
        na.rm = TRUE
      )
    }
  }
  
  ev <- eigen(med, symmetric = TRUE)
  eps <- .Machine$double.eps^(2 / 3)
  
  ev$vectors %*%
    diag(pmax(ev$values, eps), nrow = dimension, ncol = dimension) %*%
    t(ev$vectors)
}

create_het_bands <- function(mu, Sigma, fitted_model, E0, I0) {
  if (is.null(Sigma)) {
    return(NULL)
  }
  
  draws <- MASS::mvrnorm(n = N_SAMPLES, mu = mu, Sigma = Sigma)
  M <- matrix(NA_real_, nrow = length(times_eval), ncol = N_SAMPLES)
  
  for (j in seq_len(N_SAMPLES)) {
    R0_j <- exp(draws[j, 1])
    v_j  <- exp(draws[j, 2])
    t0_j <- exp(draws[j, 3])
    c2_j <- expit(draws[j, 4])
    
    comp_j <- if (fitted_model == "gamma") {
      create_gamma_components(R0_j, v_j, t0_j, c2_j)
    } else {
      create_lognormal_components(R0_j, v_j, t0_j, c2_j)
    }
    
    init_j <- initDistr(q = comp_j$q, K = K, E0 = E0, I0 = I0, N = N)
    
    M[, j] <- daily_incidence_het(
      params_vec = comp_j$params,
      init_state = init_j,
      times = times_full,
      K = K
    ) %>%
      dplyr::filter(time %in% times_eval) %>%
      dplyr::pull(Inc)
  }
  
  tibble::tibble(
    time = times_eval,
    lo = apply(M, 1, stats::quantile, 0.025, na.rm = TRUE),
    hi = apply(M, 1, stats::quantile, 0.975, na.rm = TRUE)
  )
}

create_hom_bands <- function(mu, Sigma, init_state) {
  if (is.null(Sigma)) {
    return(NULL)
  }
  
  draws <- MASS::mvrnorm(n = N_SAMPLES, mu = mu, Sigma = Sigma)
  M <- matrix(NA_real_, nrow = length(times_eval), ncol = N_SAMPLES)
  
  for (j in seq_len(N_SAMPLES)) {
    params_j <- create_hom_params(
      R0 = exp(draws[j, 1]),
      t0 = exp(draws[j, 2]),
      c2 = expit(draws[j, 3])
    )
    
    M[, j] <- daily_incidence_hom(
      params_vec = params_j,
      init_state = init_state,
      times = times_full
    ) %>%
      dplyr::filter(time %in% times_eval) %>%
      dplyr::pull(Inc)
  }
  
  tibble::tibble(
    time = times_eval,
    lo = apply(M, 1, stats::quantile, 0.025, na.rm = TRUE),
    hi = apply(M, 1, stats::quantile, 0.975, na.rm = TRUE)
  )
}

Sigma_gam <- median_covariance(gamma_covariances, 4)
Sigma_log <- median_covariance(lognormal_covariances, 4)
Sigma_hom <- median_covariance(hom_covariances, 3)

bands_gam <- create_het_bands(
  mu = c(log(fit_gam$R0), log(fit_gam$v), log(fit_gam$t0), logit(fit_gam$c2)),
  Sigma = Sigma_gam,
  fitted_model = "gamma",
  E0 = E0,
  I0 = I0
)

bands_log <- create_het_bands(
  mu = c(log(fit_log$R0), log(fit_log$v), log(fit_log$t0), logit(fit_log$c2)),
  Sigma = Sigma_log,
  fitted_model = "lognormal",
  E0 = E0,
  I0 = I0
)

if (!is.null(fit_hom)) {
  bands_hom <- create_hom_bands(
    mu = c(log(fit_hom$R0), log(fit_hom$t0), logit(fit_hom$c2)),
    Sigma = Sigma_hom,
    init_state = initial_state
  )
} else {
  bands_hom <- NULL
}

# -----------------------------------------------------------------------------
# 5.10 Prediction plot
# -----------------------------------------------------------------------------

if (TRUTH_DIST == "GAMMA") {
  label_gam <- "Gamma (correct)"
  label_log <- "Lognormal (misspecified)"
} else {
  label_gam <- "Gamma (misspecified)"
  label_log <- "Lognormal (correct)"
}

bands_gam_plot <- if (!is.null(bands_gam)) dplyr::filter(bands_gam, time > T_fit) else NULL
bands_log_plot <- if (!is.null(bands_log)) dplyr::filter(bands_log, time > T_fit) else NULL
bands_hom_plot <- if (!is.null(bands_hom)) dplyr::filter(bands_hom, time > T_fit) else NULL

ymax <- max(
  c(
    dat_fit$reports,
    traj_gam$Inc,
    traj_log$Inc,
    traj_true$Inc,
    if (!is.null(traj_hom)) traj_hom$Inc else numeric(0)
  ),
  na.rm = TRUE
)

p <- ggplot2::ggplot() +
  { if (!is.null(bands_hom_plot))
    ggplot2::geom_ribbon(
      data = bands_hom_plot,
      ggplot2::aes(x = time, ymin = lo, ymax = hi),
      fill = "forestgreen",
      alpha = 0.18,
      colour = NA
    )
  } +
  { if (!is.null(bands_log_plot))
    ggplot2::geom_ribbon(
      data = bands_log_plot,
      ggplot2::aes(x = time, ymin = lo, ymax = hi),
      fill = "darkorange",
      alpha = 0.18,
      colour = NA
    )
  } +
  { if (!is.null(bands_gam_plot))
    ggplot2::geom_ribbon(
      data = bands_gam_plot,
      ggplot2::aes(x = time, ymin = lo, ymax = hi),
      fill = "steelblue",
      alpha = 0.22,
      colour = NA
    )
  } +
  ggplot2::geom_vline(
    xintercept = T_fit,
    linetype = "dashed",
    colour = "#7B2D8B",
    linewidth = 0.6
  ) +
  ggplot2::annotate(
    "text",
    x = T_fit + 1,
    y = 0.90 * ymax,
    label = "Forecast begins",
    colour = "#7B2D8B",
    angle = 90,
    vjust = -0.4,
    hjust = 1,
    size = 3.2
  ) +
  ggplot2::geom_line(
    data = traj_gam,
    ggplot2::aes(time, Inc, colour = label_gam),
    linewidth = 1.1
  ) +
  ggplot2::geom_line(
    data = traj_log,
    ggplot2::aes(time, Inc, colour = label_log),
    linewidth = 1.1
  ) +
  { if (!is.null(traj_hom))
    ggplot2::geom_line(
      data = traj_hom,
      ggplot2::aes(time, Inc, colour = "Homogeneous"),
      linewidth = 1.1
    )
  } +
  ggplot2::geom_line(
    data = traj_true,
    ggplot2::aes(time, Inc, colour = "Truth"),
    linewidth = 1.0,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    data = dat_fit,
    ggplot2::aes(time, reports),
    colour = "black",
    size = 0.9,
    alpha = 0.60,
    shape = 16
  ) +
  ggplot2::scale_colour_manual(
    name = "Model",
    values = stats::setNames(
      c("steelblue", "darkorange", "forestgreen", "black"),
      c(label_gam, label_log, "Homogeneous", "Truth")
    )
  ) +
  ggplot2::scale_x_continuous(breaks = seq(0, tfinal_forecast, by = 50)) +
  ggplot2::labs(
    title = paste0("Prediction trajectories (Truth: ", TRUTH_DIST, ", CV = ", round(TRUTH_CV, 3), ")"),
    subtitle = paste0("Fitted to days 1-", T_fit, "; forecast to day ", tfinal_forecast),
    x = "Time (days)",
    y = "Daily cases"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey40"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.background = ggplot2::element_rect(
      fill = "white",
      colour = "grey80",
      linewidth = 0.3
    ),
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 8, face = "bold")
  )

print(p)

ggplot2::ggsave(
  filename = paste0("prediction_trajectories_with_hom_", SIMULATE_WITH, ".pdf"),
  plot = p,
  width = 10,
  height = 5.5
)
