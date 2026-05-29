# data-raw/run_sensitivity_sif_095_ml_v_cv.R
#
# SIF sensitivity: ML vs CV at fixed rho = 0.95.
#
# Both runs use: water covariate, no RSR, rho fixed at 0.95.
# Purpose: direct comparison of phi selected by ML vs CV at a rho value
# where ML previously showed agreement with CV (rho=0.95), to confirm
# the disagreement is driven by the near-intrinsic regime rather than
# a fundamental ML/CV difference.
#
# Outputs:
#   results_sif_rho095_ml  -- ML-tuned phi, rho=0.95 fixed
#   results_sif_rho095_cv  -- CV-tuned phi, rho=0.95 fixed

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- FALSE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)
message(sprintf("Using %d workers", parallel::detectCores() - 1L))

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
# 2. Fixed rho=0.95 Q functions
# ------------------------------------------------------------------------------

# Spatial-only Q at rho=0.95 -- used by tune_ml (handles augmentation internally)
Q_fun_095 <- local({
  S <- Matrix::Diagonal(nrow(W_queen)) - 0.95 * W_queen
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
})

# Augmented Q at rho=0.95 -- used by tune_cv (needs full augmented system)
Q_fun_095_aug <- local({
  S   <- Matrix::Diagonal(nrow(W_queen)) - 0.95 * W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  Q_a <- Matrix::forceSymmetric(
    Matrix::bdiag(Q, lambda_beta * Matrix::Diagonal(q))
  )
  function(theta) list(Q = Q_a, log_det_Q = NULL)
})

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_sif), 10L)

# ------------------------------------------------------------------------------
# 3. Shared fit helper (no performance metrics)
# ------------------------------------------------------------------------------

fit_cov <- function(Q_spatial, phi, rho_val, tuning_label, history = NULL) {
  Q_aug <- Matrix::forceSymmetric(
    Matrix::bdiag(Q_spatial, lambda_beta * Matrix::Diagonal(q))
  )
  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = Q_aug, phi = phi, solver = "cholesky"
  )
  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]
  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))
  list(
    run_name              = sprintf("sif_rho095_%s", tuning_label),
    tags                  = list(tuning     = tuning_label,
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "none",
                                 W          = "queen",
                                 rho        = 0.95),
    timestamp             = Sys.time(),
    history               = history,
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = phi,
    rho_opt               = 0.95,
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
# 5. ML run -- rho=0.95 fixed
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_rho095_ml")) {
  message("\n== SIF ML: water, rho=0.95 fixed ==")

  tuned_ml <- fastblm::tune_ml(
    y       = y_sif,
    A       = A,
    Q_fun   = Q_fun_095,
    X_fixed = X_obs_water,
    verbose = TRUE
  )

  message(sprintf("  ML: rho=0.95 (fixed)  phi=%.4f  sigma2e=%.4g",
                  tuned_ml$phi, tuned_ml$sigma2e))

  results_sif_rho095_ml <- fit_cov(
    Q_spatial     = tuned_ml$Q,
    phi           = tuned_ml$phi,
    rho_val       = 0.95,
    tuning_label  = "ml",
    history       = tuned_ml$history
  )

  usethis::use_data(results_sif_rho095_ml, overwrite = TRUE)
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned_ml$phi, tuned_ml$sigma2e))
}

# ------------------------------------------------------------------------------
# 6. CV run -- rho=0.95 fixed
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_rho095_cv")) {
  message("\n== SIF CV: water, rho=0.95 fixed ==")

  tuned_cv <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_095_aug,
    theta_init = numeric(0),
    k          = 10L,
    solver     = "cholesky",
    folds      = fold_assignments,
    parallel   = TRUE,
    verbose    = TRUE
  )

  message(sprintf("  CV: rho=0.95 (fixed)  phi=%.4f  sigma2e=%.4g",
                  tuned_cv$phi, tuned_cv$sigma2e))

  results_sif_rho095_cv <- fit_cov(
    Q_spatial    = tuned_cv$Q[seq_len(p), seq_len(p)],
    phi          = tuned_cv$phi,
    rho_val      = 0.95,
    tuning_label = "cv",
    history      = tuned_cv$history
  )

  usethis::use_data(results_sif_rho095_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned_cv$phi, tuned_cv$sigma2e))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF rho=0.95 ML vs CV results saved to data/.")
