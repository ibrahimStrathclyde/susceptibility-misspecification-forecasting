# =============================================================================
# TWO-EPIDEMIC MISSPECIFICATION ANALYSIS  (LA_SM version)
#
# Script summary:
#   - All Discretize_gamma_mult  calls -> Discretize_gamma_LA_SM
#   - All Discretize_lognormal_mult calls -> Discretize_lognormal_LA_SM
#   - Defaults: calibration="log_affine", rep="condexp", spacing="equal"
#   - Removed duplicate create_epi_panel / prediction blocks
#   - CSV export section included at end
# =============================================================================

library(deSolve)
library(tidyverse)
library(MASS)
library(gridExtra)
library(grid)
library(cowplot)

# --- Reproducibility ---
SEED <- 20251217
set.seed(SEED)

# Source the updated helper functions (LA_SM version)
source("R/1_helpers.R")

if (!exists("logit")) logit <- function(p) log(p / (1 - p))
if (!exists("expit")) expit <- function(x) 1 / (1 + exp(-x))

cat("=== TWO-EPIDEMIC MISSPECIFICATION ANALYSIS (LA_SM version) ===\n")

# =============================================================================
# STEP 1: FIXED PARAMETERS AND SETUP
# =============================================================================

N <- 100000
K <- 60
n_replicates <- 10

# Epidemiological parameters
R0_spec    <- 3.0
delta_spec <- 1 / 5.5
rho_spec   <- 0.5
gamma_spec <- 1 / 4

# Intervention parameters (shared)
t0_spec       <- 15
t1_spec       <- 20
t2_spec       <- 99
t3_spec       <- 100
tfinal_spec   <- 100L
c_value1_spec <- 1
c_value2_spec <- 0.3
c_value3_spec <- 1

times_full <- 0:tfinal_spec

# Initial conditions for two epidemics
E0_1 <- 1000;  I0_1 <- 400
E0_2 <- 100;    I0_2 <- 40

# Simulation distribution and CV
SIMULATE_WITH <- "gamma"     # "gamma" or "lognormal"
CV_value      <- 1.414

cat("\nSimulation parameters:\n")
cat("- Simulating with:", SIMULATE_WITH, "distribution\n")
cat("- CV:", round(CV_value, 3), "\n")
cat("- R0:", R0_spec, "\n")
cat("- NPI strength:", c_value2_spec, "\n")
cat("- Number of replicates:", n_replicates, "\n")
cat("- Two epidemics with different initial conditions\n")
cat("  - Epidemic 1: E0=", E0_1, ", I0=", I0_1, "\n", sep = "")
cat("  - Epidemic 2: E0=", E0_2, ", I0=", I0_2, "\n", sep = "")

# =============================================================================
# STEP 1b: GLOBALS REQUIRED BY THE HOMOGENEOUS FIT
#
#  fit3_hom_2epic_loglikwithNPI() reads v_spec, initial_state_1,
#  initial_state_2, and the t*_spec / c_value*_spec / rate parameters from
#  .GlobalEnv. v_spec = 0 forces the homogeneous case.
#
#  This block must run BEFORE Step 3; otherwise the homogeneous fit will fail
#  with "object 'v_spec' not found".
# =============================================================================

v_spec <- 0
assign("v_spec", v_spec, envir = .GlobalEnv)

initial_state_1 <- c(S = N - E0_1 - I0_1, E = E0_1, I = I0_1, R = 0, C = 0)
initial_state_2 <- c(S = N - E0_2 - I0_2, E = E0_2, I = I0_2, R = 0, C = 0)
assign("initial_state_1", initial_state_1, envir = .GlobalEnv)
assign("initial_state_2", initial_state_2, envir = .GlobalEnv)

for (.nm in c("N", "t0_spec", "t1_spec", "t2_spec", "t3_spec",
              "c_value1_spec", "c_value2_spec", "c_value3_spec",
              "rho_spec", "delta_spec", "gamma_spec", "tfinal_spec")) {
  assign(.nm, get(.nm), envir = .GlobalEnv)
}
rm(.nm)

# =============================================================================
# STEP 2: SIMULATE DATA FOR TWO EPIDEMICS
# =============================================================================

cat("\n=== SIMULATING DATA FOR TWO EPIDEMICS ===\n")

simulated_datasets <- list()

pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
for (i in 1:n_replicates) {
  setTxtProgressBar(pb, i)
  
  if (SIMULATE_WITH == "gamma") {
    sim_result_1 <- simulate_cases_hetsus_model_gamma(
      alpha = 1 / (CV_value^2), K = K, R0 = R0_spec, delta = delta_spec,
      rho = rho_spec, gamma = gamma_spec, N = N, E0 = E0_1, I0 = I0_1,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
    sim_result_2 <- simulate_cases_hetsus_model_gamma(
      alpha = 1 / (CV_value^2), K = K, R0 = R0_spec, delta = delta_spec,
      rho = rho_spec, gamma = gamma_spec, N = N, E0 = E0_2, I0 = I0_2,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
  } else {
    sim_result_1 <- simulate_cases_hetsus_model_lognormal(
      CV = CV_value, K = K, R0 = R0_spec, delta = delta_spec,
      rho = rho_spec, gamma = gamma_spec, N = N, E0 = E0_1, I0 = I0_1,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
    sim_result_2 <- simulate_cases_hetsus_model_lognormal(
      CV = CV_value, K = K, R0 = R0_spec, delta = delta_spec,
      rho = rho_spec, gamma = gamma_spec, N = N, E0 = E0_2, I0 = I0_2,
      t0 = t0_spec, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
      c_value1 = c_value1_spec, c_value2 = c_value2_spec, c_value3 = c_value3_spec,
      tfinal = tfinal_spec
    )
  }
  
  simulated_datasets[[i]] <- list(
    data1 = sim_result_1$sim_data,
    data2 = sim_result_2$sim_data
  )
}
close(pb)
cat("\nData simulation complete for both epidemics!\n")

# =============================================================================
# STEP 3: FIT BOTH MODELS TO TWO-EPIDEMIC DATA
# =============================================================================

cat("\n=== FITTING MODELS TO TWO-EPIDEMIC DATA ===\n")

results_gamma     <- NULL
results_lognormal <- NULL
results_homogeneous <- NULL
gamma_covariances     <- list()
lognormal_covariances <- list()
hom_covariances       <- list()

pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
for (i in 1:n_replicates) {
  setTxtProgressBar(pb, i)
  current_data1 <- simulated_datasets[[i]]$data1
  current_data2 <- simulated_datasets[[i]]$data2
  
  # ---------- FIT GAMMA ----------
  xval_gam <- Discretize_gamma_LA_SM(
    n_groups = K, alpha = 1 / (CV_value^2), beta = 1 / (CV_value^2),
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  init_het_gamma_1 <- initDistr(q = xval_gam$q, K = K, E0 = E0_1, I0 = I0_1, N = N)
  init_het_gamma_2 <- initDistr(q = xval_gam$q, K = K, E0 = E0_2, I0 = I0_2, N = N)
  
  z_mle_gamma <- tryCatch({
    fit4_seir_2epic_poisson.loglik(
      dat1 = current_data1, dat2 = current_data2,
      init_het_1 = init_het_gamma_1, init_het_2 = init_het_gamma_2
    )
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
  xval_lognorm <- Discretize_lognormal_LA_SM(
    n_groups = K, CV = CV_value,
    spacing = "equal", rep = "condexp", calibration = "log_affine",
    print_message = FALSE
  )
  init_het_lognorm_1 <- initDistr(q = xval_lognorm$q, K = K, E0 = E0_1, I0 = I0_1, N = N)
  init_het_lognorm_2 <- initDistr(q = xval_lognorm$q, K = K, E0 = E0_2, I0 = I0_2, N = N)
  
  z_mle_lognormal <- tryCatch({
    fit4_seir_2epic_poisson.loglik_lognormal(
      dat1 = current_data1, dat2 = current_data2,
      init_het_1 = init_het_lognorm_1, init_het_2 = init_het_lognorm_2
    )
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
  
  # ---------- FIT HOMOGENEOUS (joint, reference baseline) ----------
  # Joint fit over both epidemics with v fixed at 0. The objective function
  # (f_optim_reducedm.poisloglik_2epi_withNPI inside fit3_hom_2epic_...) reads
  # initial_state_1 / initial_state_2 / v_spec from .GlobalEnv (set in Step 1b).
  z_mle_hom <- tryCatch({
    fit3_hom_2epic_loglikwithNPI(dat1 = current_data1, dat2 = current_data2)
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
          
          # Homogeneous joint fit has 3 parameters: (log R0, log t0, logit c_value2)
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
cat("\nModel fitting complete for two-epidemic data!\n")

# =============================================================================
# STEP 4: ANALYZE RESULTS
# =============================================================================

cat("\n=== ANALYZING RESULTS ===\n")

valid_gamma       <- results_gamma       %>% filter(hess_pd == TRUE, convergence == 0)
valid_lognormal   <- results_lognormal   %>% filter(hess_pd == TRUE, convergence == 0)
valid_homogeneous <- if (!is.null(results_homogeneous)) {
  results_homogeneous %>% filter(hess_pd == TRUE, convergence == 0)
} else {
  results_homogeneous[0, , drop = FALSE]
}

cat("\nConvergence summary (Two Epidemics):\n")
cat("- Gamma model: ",       nrow(valid_gamma),       "out of", n_replicates, "converged (PD Hessian)\n")
cat("- Lognormal model: ",   nrow(valid_lognormal),   "out of", n_replicates, "converged (PD Hessian)\n")
cat("- Homogeneous model: ", nrow(valid_homogeneous), "out of", n_replicates, "converged (PD Hessian)\n")



# --- Export combined valid results for external plotting ---
valid_gamma$Model       <- "Gamma"
valid_lognormal$Model   <- "Lognormal"
valid_homogeneous$Model <- "Homogeneous"

combined_data_2epi <- bind_rows(valid_gamma, valid_lognormal, valid_homogeneous)
combined_data_2epi$Model <- factor(combined_data_2epi$Model, levels = c("Gamma", "Lognormal", "Homogeneous"))

write.csv(combined_data_2epi,
          paste0("combined_valid_two_epi_sim_", SIMULATE_WITH, ".csv"),
          row.names = FALSE)
cat("Saved: combined_valid_two_epi_sim_", SIMULATE_WITH, ".csv\n", sep = "")

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



# --- Two epidemics (load from the other directory first) ---
tab_gamma_2epi     <- summarise_model(valid_gamma,       "Gamma")
tab_lognormal_2epi <- summarise_model(valid_lognormal,   "Lognormal")
tab_hom_2epi       <- summarise_model(valid_homogeneous, "Homogeneous")

summary_2epi <- bind_rows(tab_gamma_2epi, tab_lognormal_2epi, tab_hom_2epi)

summary_2epi$Median_SD <- ifelse(
  is.na(summary_2epi$Median),
  "---",
  sprintf("%.3f (%.3f)", summary_2epi$Median, summary_2epi$SD)
)
summary_2epi$Width_fmt <- ifelse(
  is.na(summary_2epi$CI_Width),
  "---",
  sprintf("%.3f", summary_2epi$CI_Width)
)
summary_2epi$Cov_fmt <- ifelse(
  is.na(summary_2epi$Coverage),
  "---",
  sprintf("%.1f", summary_2epi$Coverage)
)

cat("\n=== TWO EPIDEMICS PARAMETER SUMMARY ===\n")
print(summary_2epi %>% dplyr::select(Model, Parameter, Median_SD, Width_fmt, Cov_fmt))




summary_gamma_2epi <- valid_gamma %>%
  summarise(
    Model = "Gamma", Scenario = "Two Epidemics",
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

summary_lognormal_2epi <- valid_lognormal %>%
  summarise(
    Model = "Lognormal", Scenario = "Two Epidemics",
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

summary_table_2epi <- bind_rows(summary_gamma_2epi, summary_lognormal_2epi)

# --- Homogeneous (joint) summary: v-related columns are NA ---
if (nrow(valid_homogeneous) > 0) {
  summary_homogeneous_2epi <- valid_homogeneous %>%
    summarise(
      Model = "Homogeneous", Scenario = "Two Epidemics",
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
  summary_table_2epi <- bind_rows(summary_table_2epi, summary_homogeneous_2epi)
}

cat("\nSUMMARY RESULTS (TWO EPIDEMICS):\n")
cat("Data simulated with:", SIMULATE_WITH, "distribution (CV =", CV_value, ")\n\n")
print(summary_table_2epi %>% dplyr::select(Model, R0_mean, v_mean, t0_mean, c_value2_mean))
cat("\nRELATIVE BIAS (%):\n")
print(summary_table_2epi %>% dplyr::select(Model, R0_rel_bias, v_rel_bias, c_value2_rel_bias))
cat("\nCOVERAGE:\n")
print(summary_table_2epi %>% dplyr::select(Model, R0_coverage, v_coverage, t0_coverage, c_value2_coverage))
cat("\nIDENTIFIABILITY:\n")
print(summary_table_2epi %>% dplyr::select(Model, median_condition, median_v_c_corr))

# Save results
write.csv(results_gamma,      paste0("results_gamma_2epi_sim_",     SIMULATE_WITH, ".csv"), row.names = FALSE)
write.csv(results_lognormal,  paste0("results_lognormal_2epi_sim_", SIMULATE_WITH, ".csv"), row.names = FALSE)
if (!is.null(results_homogeneous)) {
  write.csv(results_homogeneous, paste0("results_homogeneous_2epi_sim_", SIMULATE_WITH, ".csv"), row.names = FALSE)
}
write.csv(summary_table_2epi, paste0("summary_2epi_sim_",           SIMULATE_WITH, ".csv"), row.names = FALSE)
saveRDS(gamma_covariances,     paste0("gamma_covariances_2epi_sim_",     SIMULATE_WITH, ".rds"))
saveRDS(lognormal_covariances, paste0("lognormal_covariances_2epi_sim_", SIMULATE_WITH, ".rds"))
saveRDS(hom_covariances,       paste0("hom_covariances_2epi_sim_",       SIMULATE_WITH, ".rds"))

# =============================================================================
# STEP 5: VISUALIZATIONS  (density, bias, coverage, correlation)
# =============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

valid_gamma$Model     <- "Gamma"
valid_lognormal$Model <- "Lognormal"
combined_data <- bind_rows(valid_gamma, valid_lognormal)
combined_data$Model <- factor(combined_data$Model, levels = c("Gamma", "Lognormal"))

pal_2epi <- c("Gamma" = "lightblue", "Lognormal" = "lightgreen")

# --- Density plots ---
R0_plot <- ggplot(combined_data, aes(x = R0, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = R0_spec, linetype = "dashed", size = 1) +
  scale_fill_manual(values = pal_2epi) +
  theme_minimal() + theme(legend.position = "none",
                          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "R0", x = "R0", y = "Density")

CV_plot <- ggplot(combined_data, aes(x = v, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = CV_value, linetype = "dashed", size = 1) +
  scale_fill_manual(values = pal_2epi) +
  theme_minimal() + theme(legend.position = "none",
                          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "CV", x = "CV", y = "Density")

t0_plot <- ggplot(combined_data, aes(x = t0, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = t0_spec, linetype = "dashed", size = 1) +
  scale_fill_manual(values = pal_2epi) +
  theme_minimal() + theme(legend.position = "bottom",
                          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "t0", x = "t0 (days)", y = "Density")

NPI_plot <- ggplot(combined_data, aes(x = c_value2, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = c_value2_spec, linetype = "dashed", size = 1) +
  scale_fill_manual(values = pal_2epi) +
  theme_minimal() + theme(legend.position = "bottom",
                          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "NPI", x = "NPI", y = "Density")

density_combined <- grid.arrange(
  R0_plot, CV_plot, t0_plot, NPI_plot, ncol = 2, nrow = 2,
  top = textGrob(
    paste("Parameter Distributions - Two Epidemics\nSimulated with", toupper(SIMULATE_WITH)),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
ggsave(paste0("param_distributions_2epi_sim_", SIMULATE_WITH, ".png"),
       density_combined, width = 10, height = 10)

# --- Bias plot ---
bias_data <- data.frame(
  Parameter = rep(c("R0", "CV", "t0", "NPI"), 2),
  Model = rep(c("Gamma", "Lognormal"), each = 4),
  Relative_Bias = c(
    summary_gamma_2epi$R0_rel_bias, summary_gamma_2epi$v_rel_bias,
    summary_gamma_2epi$t0_bias / t0_spec * 100, summary_gamma_2epi$c_value2_rel_bias,
    summary_lognormal_2epi$R0_rel_bias, summary_lognormal_2epi$v_rel_bias,
    summary_lognormal_2epi$t0_bias / t0_spec * 100, summary_lognormal_2epi$c_value2_rel_bias
  )
)
bias_data$Model <- factor(bias_data$Model, levels = c("Gamma", "Lognormal"))

bias_plot <- ggplot(bias_data, aes(x = Parameter, y = Relative_Bias, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = pal_2epi) +
  labs(title = paste("Parameter Bias - Two Epidemics\nSimulated with", toupper(SIMULATE_WITH)),
       x = "Parameter", y = "Relative Bias (%)",
       subtitle = paste("True CV =", CV_value)) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                          plot.subtitle = element_text(hjust = 0.5))
ggsave(paste0("bias_comparison_2epi_sim_", SIMULATE_WITH, ".png"), bias_plot, width = 8, height = 6)

# --- Coverage plot ---
coverage_data <- data.frame(
  Parameter = rep(c("R0", "CV", "t0", "NPI"), 2),
  Model = rep(c("Gamma", "Lognormal"), each = 4),
  Coverage = c(
    summary_gamma_2epi$R0_coverage, summary_gamma_2epi$v_coverage,
    summary_gamma_2epi$t0_coverage, summary_gamma_2epi$c_value2_coverage,
    summary_lognormal_2epi$R0_coverage, summary_lognormal_2epi$v_coverage,
    summary_lognormal_2epi$t0_coverage, summary_lognormal_2epi$c_value2_coverage
  )
)
coverage_data$Model <- factor(coverage_data$Model, levels = c("Gamma", "Lognormal"))

coverage_plot <- ggplot(coverage_data, aes(x = Parameter, y = Coverage, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  scale_fill_manual(values = pal_2epi) + ylim(0, 1) +
  labs(title = paste("Coverage - Two Epidemics\nSimulated with", toupper(SIMULATE_WITH)),
       x = "Parameter", y = "Coverage",
       subtitle = "Red line = nominal 95%") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                          plot.subtitle = element_text(hjust = 0.5))
ggsave(paste0("coverage_comparison_2epi_sim_", SIMULATE_WITH, ".png"),
       coverage_plot, width = 8, height = 6)

# --- Correlation heatmaps ---
param_labels <- c(R0 = expression(R[0]), v = expression(nu),
                  c_value2 = expression(c), t0 = expression(t[0]))
params_all   <- names(param_labels)

create_corr_plot_full <- function(df, title_expr, params, labels) {
  prs <- params[params %in% names(df)]
  M   <- cor(df[, prs, drop = FALSE], use = "pairwise.complete.obs")
  corr_long <- as.data.frame(as.table(M)) %>%
    rename(Var1 = Var1, Var2 = Var2, Cor = Freq) %>%
    mutate(Var1 = factor(Var1, levels = prs), Var2 = factor(Var2, levels = prs))
  
  ggplot(corr_long, aes(x = Var2, y = Var1, fill = Cor)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.3f", Cor)), size = 4.5, fontface = "bold") +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         limits = c(-1, 1), name = "Correlation") +
    scale_x_discrete(labels = labels) + scale_y_discrete(labels = labels) +
    coord_equal() + labs(title = title_expr, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
}

p_gam_corr <- create_corr_plot_full(
  valid_gamma,
  expression("Gamma (two epidemics)"),
  params_all, param_labels
)
p_lgn_corr <- create_corr_plot_full(
  valid_lognormal,
  expression("Lognormal (two epidemics)"),
  params_all, param_labels
)

corr_side <- plot_grid(p_gam_corr, p_lgn_corr, ncol = 2)
ggsave(paste0("correlations_2epi_sim_", SIMULATE_WITH, ".png"),
       corr_side, width = 12, height = 6)
# =============================================================================
# STEP 6: TWO-EPIDEMIC PREDICTION TRAJECTORIES
#
# This block produces a single two-panel prediction figure containing:
#   1. Gamma fitted model
#   2. Lognormal fitted model
#   3. Homogeneous fitted model
#   4. Deterministic truth
#
# It removes the duplicated non-homogeneous-only plotting block and keeps one
# clean prediction workflow.
# =============================================================================

T_fit           <- 100L
tfinal_forecast <- 250L
times_full      <- 0:tfinal_forecast
times_eval      <- (T_fit + 1):tfinal_forecast
N_SAMPLES       <- 300L

TRUTH_DIST <- toupper(SIMULATE_WITH)
TRUTH_CV   <- CV_value
sw         <- SIMULATE_WITH

# -----------------------------------------------------------------------------
# 6.1 Mean observed training data across simulation replicates
# -----------------------------------------------------------------------------

dat_fit_1 <- dplyr::bind_rows(lapply(simulated_datasets, `[[`, "data1")) %>%
  dplyr::filter(time >= 1, time <= T_fit) %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(reports = mean(reports, na.rm = TRUE), .groups = "drop")

dat_fit_2 <- dplyr::bind_rows(lapply(simulated_datasets, `[[`, "data2")) %>%
  dplyr::filter(time >= 1, time <= T_fit) %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(reports = mean(reports, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# 6.2 Daily-incidence helpers
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
# 6.3 Discretisation and parameter-vector builders
# -----------------------------------------------------------------------------

disc_gamma <- function(v) {
  alpha <- 1 / (v^2)
  Discretize_gamma_LA_SM(
    n_groups = K,
    alpha = alpha,
    beta = alpha,
    spacing = "equal",
    rep = "condexp",
    calibration = "log_affine",
    print_message = FALSE
  )
}

disc_lognormal <- function(v) {
  Discretize_lognormal_LA_SM(
    n_groups = K,
    CV = v,
    spacing = "equal",
    rep = "condexp",
    calibration = "log_affine",
    print_message = FALSE
  )
}

create_gamma_components <- function(R0, v, t0, c2) {
  d <- disc_gamma(v)
  
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
  d <- disc_lognormal(v)
  
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
# 6.4 Median parameter estimates
# -----------------------------------------------------------------------------

clip01 <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)

median_fit_het <- function(df, model_label) {
  if (is.null(df) || nrow(df) == 0) {
    stop(paste("No valid", model_label, "fits are available."))
  }
  
  list(
    R0 = max(stats::median(df$R0, na.rm = TRUE), 1e-8),
    v  = max(stats::median(df$v, na.rm = TRUE),  1e-8),
    t0 = max(stats::median(df$t0, na.rm = TRUE), 1e-8),
    c2 = clip01(stats::median(df$c_value2, na.rm = TRUE))
  )
}

median_fit_hom <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  
  list(
    R0 = max(stats::median(df$R0, na.rm = TRUE), 1e-8),
    t0 = max(stats::median(df$t0, na.rm = TRUE), 1e-8),
    c2 = clip01(stats::median(df$c_value2, na.rm = TRUE))
  )
}

fit_gam <- median_fit_het(valid_gamma, "Gamma")
fit_log <- median_fit_het(valid_lognormal, "Lognormal")
fit_hom <- median_fit_hom(valid_homogeneous)

# -----------------------------------------------------------------------------
# 6.5 Deterministic truth trajectories
# -----------------------------------------------------------------------------

truth_comp <- if (SIMULATE_WITH == "gamma") {
  create_gamma_components(R0_spec, CV_value, t0_spec, c_value2_spec)
} else {
  create_lognormal_components(R0_spec, CV_value, t0_spec, c_value2_spec)
}

init_truth_1 <- initDistr(q = truth_comp$q, K = K, E0 = E0_1, I0 = I0_1, N = N)
init_truth_2 <- initDistr(q = truth_comp$q, K = K, E0 = E0_2, I0 = I0_2, N = N)

traj_true_1 <- daily_incidence_het(truth_comp$params, init_truth_1, times_full, K) %>%
  dplyr::filter(time >= 1) %>%
  dplyr::rename(Inc_truth = Inc)

traj_true_2 <- daily_incidence_het(truth_comp$params, init_truth_2, times_full, K) %>%
  dplyr::filter(time >= 1) %>%
  dplyr::rename(Inc_truth = Inc)

# -----------------------------------------------------------------------------
# 6.6 Mean fitted trajectories: Gamma, Lognormal, and Homogeneous
# -----------------------------------------------------------------------------

gam_comp <- create_gamma_components(fit_gam$R0, fit_gam$v, fit_gam$t0, fit_gam$c2)
log_comp <- create_lognormal_components(fit_log$R0, fit_log$v, fit_log$t0, fit_log$c2)

init_gam_1 <- initDistr(q = gam_comp$q, K = K, E0 = E0_1, I0 = I0_1, N = N)
init_gam_2 <- initDistr(q = gam_comp$q, K = K, E0 = E0_2, I0 = I0_2, N = N)

init_log_1 <- initDistr(q = log_comp$q, K = K, E0 = E0_1, I0 = I0_1, N = N)
init_log_2 <- initDistr(q = log_comp$q, K = K, E0 = E0_2, I0 = I0_2, N = N)

traj_gam_1 <- daily_incidence_het(gam_comp$params, init_gam_1, times_full, K) %>%
  dplyr::filter(time >= 1)

traj_gam_2 <- daily_incidence_het(gam_comp$params, init_gam_2, times_full, K) %>%
  dplyr::filter(time >= 1)

traj_log_1 <- daily_incidence_het(log_comp$params, init_log_1, times_full, K) %>%
  dplyr::filter(time >= 1)

traj_log_2 <- daily_incidence_het(log_comp$params, init_log_2, times_full, K) %>%
  dplyr::filter(time >= 1)

if (!is.null(fit_hom)) {
  params_hom <- create_hom_params(fit_hom$R0, fit_hom$t0, fit_hom$c2)
  
  traj_hom_1 <- daily_incidence_hom(params_hom, initial_state_1, times_full) %>%
    dplyr::filter(time >= 1)
  
  traj_hom_2 <- daily_incidence_hom(params_hom, initial_state_2, times_full) %>%
    dplyr::filter(time >= 1)
} else {
  traj_hom_1 <- NULL
  traj_hom_2 <- NULL
}

# -----------------------------------------------------------------------------
# 6.7 Forecast bands
# -----------------------------------------------------------------------------

median_covariance <- function(cov_list, dimension) {
  mats <- Filter(
    function(M) is.matrix(M) &&
      all(dim(M) == c(dimension, dimension)) &&
      all(is.finite(M)),
    cov_list
  )
  
  if (!length(mats)) return(NULL)
  
  med <- matrix(NA_real_, dimension, dimension)
  for (i in seq_len(dimension)) {
    for (j in seq_len(dimension)) {
      med[i, j] <- stats::median(
        vapply(mats, function(M) M[i, j], numeric(1)),
        na.rm = TRUE
      )
    }
  }
  
  ev  <- eigen(med, symmetric = TRUE)
  eps <- .Machine$double.eps^(2 / 3)
  
  ev$vectors %*%
    diag(pmax(ev$values, eps), nrow = dimension, ncol = dimension) %*%
    t(ev$vectors)
}

create_bands_het <- function(mu, Sigma, dist, E0, I0) {
  if (is.null(Sigma)) return(NULL)
  
  draws <- MASS::mvrnorm(n = N_SAMPLES, mu = mu, Sigma = Sigma)
  M <- matrix(NA_real_, nrow = length(times_eval), ncol = N_SAMPLES)
  
  for (j in seq_len(N_SAMPLES)) {
    R0_j <- exp(draws[j, 1])
    v_j  <- exp(draws[j, 2])
    t0_j <- exp(draws[j, 3])
    c2_j <- expit(draws[j, 4])
    
    comp_j <- if (dist == "gamma") {
      create_gamma_components(R0_j, v_j, t0_j, c2_j)
    } else {
      create_lognormal_components(R0_j, v_j, t0_j, c2_j)
    }
    
    init_j <- initDistr(q = comp_j$q, K = K, E0 = E0, I0 = I0, N = N)
    
    M[, j] <- daily_incidence_het(comp_j$params, init_j, times_full, K) %>%
      dplyr::filter(time %in% times_eval) %>%
      dplyr::pull(Inc)
  }
  
  tibble::tibble(
    time = times_eval,
    lo = apply(M, 1, stats::quantile, 0.025, na.rm = TRUE),
    hi = apply(M, 1, stats::quantile, 0.975, na.rm = TRUE)
  )
}

create_bands_hom <- function(mu, Sigma, init_state) {
  if (is.null(Sigma)) return(NULL)
  
  draws <- MASS::mvrnorm(n = N_SAMPLES, mu = mu, Sigma = Sigma)
  M <- matrix(NA_real_, nrow = length(times_eval), ncol = N_SAMPLES)
  
  for (j in seq_len(N_SAMPLES)) {
    R0_j <- exp(draws[j, 1])
    t0_j <- exp(draws[j, 2])
    c2_j <- expit(draws[j, 3])
    
    params_j <- create_hom_params(R0_j, t0_j, c2_j)
    
    M[, j] <- daily_incidence_hom(params_j, init_state, times_full) %>%
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

mu_gam <- c(log(fit_gam$R0), log(fit_gam$v), log(fit_gam$t0), logit(fit_gam$c2))
mu_log <- c(log(fit_log$R0), log(fit_log$v), log(fit_log$t0), logit(fit_log$c2))

bands_gam_1 <- create_bands_het(mu_gam, Sigma_gam, "gamma",     E0_1, I0_1)
bands_gam_2 <- create_bands_het(mu_gam, Sigma_gam, "gamma",     E0_2, I0_2)
bands_log_1 <- create_bands_het(mu_log, Sigma_log, "lognormal", E0_1, I0_1)
bands_log_2 <- create_bands_het(mu_log, Sigma_log, "lognormal", E0_2, I0_2)

if (!is.null(fit_hom)) {
  mu_hom <- c(log(fit_hom$R0), log(fit_hom$t0), logit(fit_hom$c2))
  
  bands_hom_1 <- create_bands_hom(mu_hom, Sigma_hom, initial_state_1)
  bands_hom_2 <- create_bands_hom(mu_hom, Sigma_hom, initial_state_2)
} else {
  bands_hom_1 <- NULL
  bands_hom_2 <- NULL
}

# -----------------------------------------------------------------------------
# 6.8 Plot labels
# -----------------------------------------------------------------------------

if (TRUTH_DIST == "GAMMA") {
  label_gam <- "Gamma (correct)"
  label_log <- "Lognormal (misspecified)"
} else {
  label_gam <- "Gamma (misspecified)"
  label_log <- "Lognormal (correct)"
}

# -----------------------------------------------------------------------------
# 6.9 Panel builder
# -----------------------------------------------------------------------------

create_prediction_panel <- function(traj_gam, traj_log, traj_hom, traj_true,
                                    bands_gam, bands_log, bands_hom,
                                    dat_fit, epi_label, show_legend = FALSE) {
  
  bg_plot <- if (!is.null(bands_gam)) dplyr::filter(bands_gam, time > T_fit) else NULL
  bl_plot <- if (!is.null(bands_log)) dplyr::filter(bands_log, time > T_fit) else NULL
  bh_plot <- if (!is.null(bands_hom)) dplyr::filter(bands_hom, time > T_fit) else NULL
  
  ymax <- max(
    c(
      dat_fit$reports,
      traj_gam$Inc,
      traj_log$Inc,
      traj_true$Inc_truth,
      if (!is.null(traj_hom)) traj_hom$Inc else numeric(0)
    ),
    na.rm = TRUE
  )
  
  ggplot2::ggplot() +
    { if (!is.null(bh_plot))
      ggplot2::geom_ribbon(
        data = bh_plot,
        ggplot2::aes(x = time, ymin = lo, ymax = hi),
        fill = "forestgreen",
        alpha = 0.15,
        colour = NA
      )
    } +
    { if (!is.null(bl_plot))
      ggplot2::geom_ribbon(
        data = bl_plot,
        ggplot2::aes(x = time, ymin = lo, ymax = hi),
        fill = "darkorange",
        alpha = 0.18,
        colour = NA
      )
    } +
    { if (!is.null(bg_plot))
      ggplot2::geom_ribbon(
        data = bg_plot,
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
      data = traj_log,
      ggplot2::aes(time, Inc, colour = label_log),
      linewidth = 1.1
    ) +
    ggplot2::geom_line(
      data = traj_gam,
      ggplot2::aes(time, Inc, colour = label_gam),
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
      ggplot2::aes(time, Inc_truth, colour = "Truth"),
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
      title = paste0(epi_label, ": prediction trajectories"),
      subtitle = paste0(
        "Truth: ", TRUTH_DIST,
        ", CV = ", round(TRUTH_CV, 3),
        "; fitted to days 1-", T_fit,
        "; forecast to day ", tfinal_forecast
      ),
      x = "Time (days)",
      y = "Daily cases"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.background = ggplot2::element_rect(
        fill = "white",
        colour = "grey80",
        linewidth = 0.3
      ),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8, face = "bold")
    )
}

# -----------------------------------------------------------------------------
# 6.10 Build and save final two-panel prediction figure
# -----------------------------------------------------------------------------

p_epi1 <- create_prediction_panel(
  traj_gam = traj_gam_1,
  traj_log = traj_log_1,
  traj_hom = traj_hom_1,
  traj_true = traj_true_1,
  bands_gam = bands_gam_1,
  bands_log = bands_log_1,
  bands_hom = bands_hom_1,
  dat_fit = dat_fit_1,
  epi_label = "Epidemic 1",
  show_legend = FALSE
)

p_epi2 <- create_prediction_panel(
  traj_gam = traj_gam_2,
  traj_log = traj_log_2,
  traj_hom = traj_hom_2,
  traj_true = traj_true_2,
  bands_gam = bands_gam_2,
  bands_log = bands_log_2,
  bands_hom = bands_hom_2,
  dat_fit = dat_fit_2,
  epi_label = "Epidemic 2",
  show_legend = TRUE
)

combined_plot <- gridExtra::arrangeGrob(p_epi1, p_epi2, ncol = 2)
print(combined_plot)
ggplot2::ggsave(
  filename = paste0("Prediction_trajectories_2epi_sim_", sw, ".png"),
  plot = combined_plot,
  width = 32,
  height = 14,
  units = "cm",
  dpi = 300
)


combined_plot <- gridExtra::arrangeGrob(p_epi1, p_epi2, ncol = 2)

grid::grid.newpage()
grid::grid.draw(combined_plot)

cat("\nTwo-epidemic prediction plot saved: ",
    paste0("Prediction_trajectories_2epi_sim_", sw, ".png"), "\n", sep = "")

# =============================================================================
# STEP 7: EXPORT TRAJECTORIES AND FORECAST BANDS
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1 Mean trajectories
# -----------------------------------------------------------------------------

write.csv(
  traj_gam_1 %>% dplyr::rename(Inc_gamma = Inc),
  paste0("2epi_traj_gamma_epi1_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  traj_gam_2 %>% dplyr::rename(Inc_gamma = Inc),
  paste0("2epi_traj_gamma_epi2_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  traj_log_1 %>% dplyr::rename(Inc_lognormal = Inc),
  paste0("2epi_traj_lognormal_epi1_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  traj_log_2 %>% dplyr::rename(Inc_lognormal = Inc),
  paste0("2epi_traj_lognormal_epi2_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  traj_true_1,
  paste0("2epi_traj_truth_epi1_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  traj_true_2,
  paste0("2epi_traj_truth_epi2_sim_", sw, ".csv"),
  row.names = FALSE
)

if (!is.null(traj_hom_1)) {
  write.csv(
    traj_hom_1 %>% dplyr::rename(Inc_homogeneous = Inc),
    paste0("2epi_traj_homogeneous_epi1_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(traj_hom_2)) {
  write.csv(
    traj_hom_2 %>% dplyr::rename(Inc_homogeneous = Inc),
    paste0("2epi_traj_homogeneous_epi2_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 7.2 Forecast bands
# -----------------------------------------------------------------------------

if (!is.null(bands_gam_1)) {
  write.csv(
    bands_gam_1,
    paste0("2epi_bands_gamma_epi1_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(bands_gam_2)) {
  write.csv(
    bands_gam_2,
    paste0("2epi_bands_gamma_epi2_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(bands_log_1)) {
  write.csv(
    bands_log_1,
    paste0("2epi_bands_lognormal_epi1_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(bands_log_2)) {
  write.csv(
    bands_log_2,
    paste0("2epi_bands_lognormal_epi2_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(bands_hom_1)) {
  write.csv(
    bands_hom_1,
    paste0("2epi_bands_homogeneous_epi1_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

if (!is.null(bands_hom_2)) {
  write.csv(
    bands_hom_2,
    paste0("2epi_bands_homogeneous_epi2_sim_", sw, ".csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 7.3 Observed training data
# -----------------------------------------------------------------------------

write.csv(
  dat_fit_1,
  paste0("2epi_obs_epi1_sim_", sw, ".csv"),
  row.names = FALSE
)

write.csv(
  dat_fit_2,
  paste0("2epi_obs_epi2_sim_", sw, ".csv"),
  row.names = FALSE
)

cat("\nTwo-epidemic prediction exports complete.\n")