# data-raw/run_sif_water_rho_cv.R
#
# SIF downscaling: water covariate, rho CV-tuned, no RSR.
# Intended for direct comparison with results_sif_water_ml (same model,
# ML tuning instead of CV) and results_sif_rsr (same tuning, adds RSR).
#
# Outputs saved to data/ via usethis::use_data():
#   results_sif_water_rho_cv

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
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_sif <- goebel2026::soundings_augmented$SIF_757nm

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

make_Q_fun_rho_aug <- function(W, q, lambda) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q   <- Matrix::drop0(Q)
    Q_a <- Matrix::forceSymmetric(
      Matrix::bdiag(Q, lambda * Matrix::Diagonal(q))
    )
    list(Q = Q_a, log_det_Q = NULL)
  }
}

Q_fun_rho_aug <- make_Q_fun_rho_aug(W_queen, q, lambda_beta)

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_sif), 10L)

# ------------------------------------------------------------------------------
# 3. Skip helper
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
# 4. Run
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_water_rho_cv")) {
  message("\n== SIF: water covariate, rho CV-tuned, no RSR ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    k          = 10L,
    solver     = "cholesky",
    folds      = fold_assignments,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  # Final fit
  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)

  A_pred  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se      <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))

  results_sif_water_rho_cv <- list(
    run_name              = "sif_water_rho_cv",
    tags                  = list(tuning     = "cv",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "none",
                                 W          = "queen",
                                 rho        = rho_hat),
    timestamp             = Sys.time(),
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    cv_curve              = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )

  usethis::use_data(results_sif_water_rho_cv, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, fit$sigma2e))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF water rho CV results saved to data/.")
