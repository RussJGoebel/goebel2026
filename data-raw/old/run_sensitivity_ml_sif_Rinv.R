# data-raw/run_sensitivity_ml_sif_rinv.R
#
# Supplement: ML tuning for SIF with per-sounding inverse noise covariance.
#
# Mirrors results_sif_water_ml but passes R_inv (per-sounding OCO-2
# uncertainty) to tune_ml. The goal is to check whether incorporating
# heteroskedastic noise brings ML phi selection in line with CV, which
# would suggest the ML/CV discrepancy on SIF is due to noise
# misspecification rather than anything about the spatial structure.
#
# Runs:
#   1. Water covariate, rho ML-tuned,    R_inv supplied
#   2. Water covariate, rho=0.99 fixed,  R_inv supplied
#   3. No covariates,   rho ML-tuned,    R_inv supplied
#
# Outputs saved to data/ via usethis::use_data():
#   results_sif_water_ml_rinv
#   results_sif_water_ml_rho1_rinv
#   results_sif_nocov_ml_rinv

library(fastblm)
library(goebel2026)
library(Matrix)
library(usethis)

FORCE_RERUN <- FALSE

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_sif    <- goebel2026::setup_sif

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_sif  <- goebel2026::soundings_augmented$SIF_757nm
R_inv  <- d_sif$R_inv

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

message(sprintf("n=%d  p=%d  q=%d", length(y_sif), p, q))
message(sprintf("R_inv class: %s", class(R_inv)))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

make_Q_fun_rho <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q   <- Matrix::drop0(Q)
    list(Q = Q, log_det_Q = NULL)
  }
}

make_Q_fun_fixed <- function(W, rho = 0.99) {
  S <- Matrix::Diagonal(nrow(W)) - rho * W
  Q <- Matrix::forceSymmetric(Matrix::crossprod(S))
  Q <- Matrix::drop0(Q)
  function(theta) list(Q = Q, log_det_Q = NULL)
}

Q_fun_rho   <- make_Q_fun_rho(W_queen)
Q_fun_fixed <- make_Q_fun_fixed(W_queen, rho = 0.99)

# ------------------------------------------------------------------------------
# 3. Shared fit helpers
# ------------------------------------------------------------------------------

fit_and_save_cov <- function(tuned, rho_val) {
  A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")
  Q_aug <- Matrix::forceSymmetric(
    Matrix::bdiag(tuned$Q, lambda_beta * Matrix::Diagonal(q))
  )

  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = Q_aug,
    phi    = tuned$phi,
    R_inv  = R_inv,
    solver = "cholesky"
  )

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)

  A_pred  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se      <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))

  list(
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    ml_history            = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

fit_and_save_nocov <- function(tuned, rho_val) {
  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A,
    Q      = tuned$Q,
    phi    = tuned$phi,
    R_inv  = R_inv,
    solver = "cholesky"
  )

  mu <- fit$posterior_mean
  se <- fastblm::posterior_se(fit, n_probes = 200L)

  list(
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = NULL,
    se_beta               = NULL,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    ml_history            = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

# ------------------------------------------------------------------------------
# 4. Skip helper
# ------------------------------------------------------------------------------

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name))
    return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error   = function(e) NULL,
    warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}

# ------------------------------------------------------------------------------
# 5. Runs
# ------------------------------------------------------------------------------

# --- Water covariate, rho ML-tuned, R_inv ------------------------------------

if (!.should_skip("results_sif_water_ml_rinv")) {
  message("\n== SIF ML + R_inv: water covariate, rho ML-tuned ==")

  tuned <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = X_obs_water,
    R_inv      = R_inv,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  ML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  out <- fit_and_save_cov(tuned, rho_hat)

  results_sif_water_ml_rinv <- c(
    list(run_name  = "sif_water_ml_rinv",
         tags      = list(tuning     = "ml",
                          response   = "SIF_757nm",
                          covariates = "water",
                          constraint = "none",
                          W          = "queen",
                          rho        = rho_hat,
                          R_inv      = "sif_uncertainty"),
         timestamp = Sys.time()),
    out
  )

  usethis::use_data(results_sif_water_ml_rinv, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))
}

# --- Water covariate, rho=0.99 fixed, R_inv ----------------------------------

if (!.should_skip("results_sif_water_ml_rho1_rinv")) {
  message("\n== SIF ML + R_inv: water covariate, rho=0.99 fixed ==")

  tuned_rho1 <- fastblm::tune_ml(
    y       = y_sif,
    A       = A,
    Q_fun   = Q_fun_fixed,
    X_fixed = X_obs_water,
    R_inv   = R_inv,
    verbose = TRUE
  )

  message(sprintf("  ML optimum: rho=0.99 (fixed)  phi=%.4f  sigma2e=%.4g",
                  tuned_rho1$phi, tuned_rho1$sigma2e))

  out1 <- fit_and_save_cov(tuned_rho1, rho_val = 0.99)

  results_sif_water_ml_rho1_rinv <- c(
    list(run_name  = "sif_water_ml_rho1_rinv",
         tags      = list(tuning     = "ml",
                          response   = "SIF_757nm",
                          covariates = "water",
                          constraint = "none",
                          W          = "queen",
                          rho        = 0.99,
                          R_inv      = "sif_uncertainty"),
         timestamp = Sys.time()),
    out1
  )

  usethis::use_data(results_sif_water_ml_rho1_rinv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned_rho1$phi, tuned_rho1$sigma2e))
}

# --- No covariates, rho ML-tuned, R_inv --------------------------------------

if (!.should_skip("results_sif_nocov_ml_rinv")) {
  message("\n== SIF ML + R_inv: no covariates, rho ML-tuned ==")

  tuned_nocov <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = NULL,
    R_inv      = R_inv,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )

  rho_hat_nocov <- tuned_nocov$theta[["rho"]]
  message(sprintf("  ML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_nocov, tuned_nocov$phi, tuned_nocov$sigma2e))

  out_nocov <- fit_and_save_nocov(tuned_nocov, rho_hat_nocov)

  results_sif_nocov_ml_rinv <- c(
    list(run_name  = "sif_nocov_ml_rinv",
         tags      = list(tuning     = "ml",
                          response   = "SIF_757nm",
                          covariates = "none",
                          constraint = "none",
                          W          = "queen",
                          rho        = rho_hat_nocov,
                          R_inv      = "sif_uncertainty"),
         timestamp = Sys.time()),
    out_nocov
  )

  usethis::use_data(results_sif_nocov_ml_rinv, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_nocov, tuned_nocov$phi, tuned_nocov$sigma2e))
}

# ------------------------------------------------------------------------------

message("\nSIF ML + R_inv sensitivity results saved to data/.")
