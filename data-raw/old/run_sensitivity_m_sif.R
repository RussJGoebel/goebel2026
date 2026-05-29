# data-raw/run_sensitivity_ml_sif.R
#
# Marginal likelihood tuning sensitivity -- SIF data.
#
# No ground truth; purpose is to compare phi/rho selected by ML vs CV.
# All runs: queen adjacency, no RSR, water covariate unless noted.
#
# Outputs:
#   results_sif_water_ml         -- water covariate, rho ML-tuned
#   results_sif_water_ml_rho1    -- water covariate, rho=0.99 fixed
#   results_sif_nocov_ml         -- no covariates, rho ML-tuned
#   results_sif_water_ml_zerolog -- water covariate, rho ML-tuned, logdet(Q)=0

library(fastblm)
library(goebel2026)
library(Matrix)
library(usethis)

FORCE_RERUN <- FALSE

# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared

A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water
fine_grid_buffered <- d_shared$fine_grid_buffered

y_sif <- goebel2026::soundings_augmented$SIF_757nm

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Q functions
# ------------------------------------------------------------------------------

Q_fun_rho <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  list(Q = Q, log_det_Q = NULL)
}

Q_fun_fixed_99 <- local({
  S <- Matrix::Diagonal(nrow(W_queen)) - 0.99 * W_queen
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
})

Q_fun_zerolog <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  list(Q = Q, log_det_Q = 0)
}

# ------------------------------------------------------------------------------
# 3. Fit helpers (no performance metrics)
# ------------------------------------------------------------------------------

fit_cov <- function(tuned, rho_val) {
  A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")
  Q_aug <- Matrix::forceSymmetric(
    Matrix::bdiag(tuned$Q, lambda_beta * Matrix::Diagonal(q))
  )
  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = Q_aug, phi = tuned$phi, solver = "cholesky"
  )
  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]
  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))
  list(
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

fit_nocov <- function(tuned, rho_val) {
  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A, Q = tuned$Q, phi = tuned$phi, solver = "cholesky"
  )
  mu <- fit$posterior_mean
  se <- fastblm::posterior_se(fit, n_probes = 200L)
  list(
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

# ------------------------------------------------------------------------------
# 4. Skip helper
# ------------------------------------------------------------------------------

.should_skip <- function(nm) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(nm, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s", nm)); return(TRUE)
  }
  pkg <- tryCatch(
    utils::data(list = nm, package = "goebel2026", envir = new.env()),
    warning = function(w) NULL, error = function(e) NULL
  )
  if (!is.null(pkg)) {
    message(sprintf("  skipping %s", nm)); return(TRUE)
  }
  FALSE
}

# ------------------------------------------------------------------------------
# 5. Runs
# ------------------------------------------------------------------------------

# --- 1. Water, rho ML-tuned ---------------------------------------------------

if (!.should_skip("results_sif_water_ml")) {
  message("\n== SIF ML: water, rho tuned ==")
  t <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = X_obs_water,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )
  rho_hat <- t$theta[["rho"]]
  message(sprintf("  rho=%.4f  phi=%.4f  sigma2e=%.4g", rho_hat, t$phi, t$sigma2e))
  out <- fit_cov(t, rho_hat)
  results_sif_water_ml <- c(
    list(run_name  = "sif_water_ml",
         tags      = list(tuning="ml", response="SIF_757nm", covariates="water",
                          constraint="none", W="queen", rho=rho_hat),
         timestamp  = Sys.time(),
         ml_history = t$history),
    out)
  usethis::use_data(results_sif_water_ml, overwrite = TRUE)
}

# --- 2. Water, rho=0.99 fixed -------------------------------------------------

if (!.should_skip("results_sif_water_ml_rho1")) {
  message("\n== SIF ML: water, rho=0.99 fixed ==")
  t <- fastblm::tune_ml(
    y       = y_sif,
    A       = A,
    Q_fun   = Q_fun_fixed_99,
    X_fixed = X_obs_water,
    verbose = TRUE
  )
  message(sprintf("  rho=0.99 (fixed)  phi=%.4f  sigma2e=%.4g", t$phi, t$sigma2e))
  out <- fit_cov(t, rho_val = 0.99)
  results_sif_water_ml_rho1 <- c(
    list(run_name  = "sif_water_ml_rho1",
         tags      = list(tuning="ml", response="SIF_757nm", covariates="water",
                          constraint="none", W="queen", rho=0.99),
         timestamp  = Sys.time(),
         ml_history = t$history),
    out)
  usethis::use_data(results_sif_water_ml_rho1, overwrite = TRUE)
}

# --- 3. No covariates, rho ML-tuned -------------------------------------------

if (!.should_skip("results_sif_nocov_ml")) {
  message("\n== SIF ML: no covariates, rho tuned ==")
  t <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = NULL,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )
  rho_hat <- t$theta[["rho"]]
  message(sprintf("  rho=%.4f  phi=%.4f  sigma2e=%.4g", rho_hat, t$phi, t$sigma2e))
  out <- fit_nocov(t, rho_hat)
  results_sif_nocov_ml <- c(
    list(run_name  = "sif_nocov_ml",
         tags      = list(tuning="ml", response="SIF_757nm", covariates="none",
                          constraint="none", W="queen", rho=rho_hat),
         timestamp  = Sys.time(),
         ml_history = t$history),
    out)
  usethis::use_data(results_sif_nocov_ml, overwrite = TRUE)
}

# --- 4. Water, rho ML-tuned, logdet(Q) = 0 ------------------------------------
# Tests whether the near-singularity of Q as rho -> 1 distorts phi selection.
# Setting log_det_Q = 0 drops the logdet(Q) term from the ML objective entirely,
# matching the improper-prior convention.

if (!.should_skip("results_sif_water_ml_zerolog")) {
  message("\n== SIF ML: water, rho tuned, logdet(Q)=0 ==")
  t <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_zerolog,
    X_fixed    = X_obs_water,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )
  rho_hat <- t$theta[["rho"]]
  message(sprintf("  rho=%.4f  phi=%.4f  sigma2e=%.4g", rho_hat, t$phi, t$sigma2e))
  out <- fit_cov(t, rho_hat)
  results_sif_water_ml_zerolog <- c(
    list(run_name  = "sif_water_ml_zerolog",
         tags      = list(tuning="ml_zerolog", response="SIF_757nm", covariates="water",
                          constraint="none", W="queen", rho=rho_hat),
         timestamp  = Sys.time(),
         ml_history = t$history),
    out)
  usethis::use_data(results_sif_water_ml_zerolog, overwrite = TRUE)
}

# ------------------------------------------------------------------------------

message("\nSIF ML sensitivity results saved to data/.")
