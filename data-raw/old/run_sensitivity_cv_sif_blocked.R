# data-raw/run_sensitivity_sif_cv_blocked.R
#
# SIF sensitivity: spatially blocked vs random CV, water + RSR, rho=1.
#
# Mirrors the canonical SIF run (results_sif_canonical) but uses spatially
# blocked folds. Purpose: check whether spatial blocking changes phi
# selection on real SIF data with RSR, analogous to the albedo blocked CV
# experiment (results_cv_blocked).
#
# Outputs:
#   results_sif_cv_blocked          -- blocked CV, water + RSR, rho=1
#   results_sif_cv_blocked_rho_cv   -- blocked CV, water + RSR, rho CV-tuned

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

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

# Intrinsic SAR prior, augmented
Q_fun_fixed_aug <- local({
  S   <- Matrix::Diagonal(nrow(W_queen)) - W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  Q_a <- Matrix::forceSymmetric(
    Matrix::bdiag(Q, lambda_beta * Matrix::Diagonal(q))
  )
  function(theta) list(Q = Q_a, log_det_Q = NULL)
})

# RSR constraint -- per fold
rsr_constraint <- local({
  force(X_obs_water); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs_water[train_idx, , drop = FALSE]) %*% A_train)
    C_zeros   <- matrix(0, nrow = q, ncol = q)
    cbind(C_spatial, C_zeros)
  }
})

# Full-data RSR constraint for final fit
C_aug_full <- cbind(
  as.matrix(t(X_obs_water) %*% A),
  matrix(0, nrow = q, ncol = q)
)

# Blocked folds -- try setup_sif first, fall back to setup_albedo
blocked_folds <- goebel2026::setup_sif$blocked_folds %||%
  goebel2026::setup_albedo$blocked_folds

if (is.null(blocked_folds)) stop("blocked_folds not found in setup_sif or setup_albedo")
message(sprintf("Blocked folds: %d blocks, sizes: %s",
                max(blocked_folds),
                paste(table(blocked_folds), collapse = " ")))

# ------------------------------------------------------------------------------
# 3. Skip helper
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
# 4. Run
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_cv_blocked")) {
  message("\n== SIF CV blocked: water + RSR, rho=1 ==")

  tuned <- fastblm::tune_cv(
    y             = y_sif,
    A             = A_aug,
    Q_fun         = Q_fun_fixed_aug,
    theta_init    = numeric(0),
    k             = 10L,
    solver        = "cholesky",
    folds         = blocked_folds,
    constraint    = rsr_constraint,
    log_phi_lower = log(0.001),
    parallel      = TRUE,
    verbose       = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g",
                  tuned$phi, tuned$sigma2e))

  # Final fit with RSR applied
  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))
  message(sprintf("  phi comparison -- canonical: %.4f  blocked: %.4f",
                  goebel2026::results_sif_canonical$phi, tuned$phi))

  results_sif_cv_blocked <- list(
    run_name              = "sif_cv_blocked",
    tags                  = list(tuning     = "cv_blocked",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "RSR",
                                 W          = "queen",
                                 rho        = 1,
                                 n_blocks   = max(blocked_folds)),
    timestamp             = Sys.time(),
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = 1,
    cv_curve              = tuned$history,
    fold_assignments      = blocked_folds,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )

  usethis::use_data(results_sif_cv_blocked, overwrite = TRUE)
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned$phi, fit$sigma2e))
}

# ------------------------------------------------------------------------------
# 5. Blocked CV run -- rho tuned
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_cv_blocked_rho_cv")) {
  message("\n== SIF CV blocked: water + RSR, rho CV-tuned ==")

  # Rho-tunable augmented Q fun
  Q_fun_rho_aug <- function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
    Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
    Q_a <- Matrix::forceSymmetric(
      Matrix::bdiag(Q, lambda_beta * Matrix::Diagonal(q))
    )
    list(Q = Q_a, log_det_Q = NULL)
  }

  tuned_rho <- fastblm::tune_cv(
    y             = y_sif,
    A             = A_aug,
    Q_fun         = Q_fun_rho_aug,
    theta_init    = c(rho = 0.9),
    lower         = c(rho = 0.5),
    upper         = c(rho = 0.999),
    k             = 10L,
    solver        = "cholesky",
    folds         = blocked_folds,
    constraint    = rsr_constraint,
    log_phi_lower = log(0.001),
    parallel      = TRUE,
    verbose       = TRUE
  )

  rho_hat <- tuned_rho$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned_rho$phi, tuned_rho$sigma2e))

  fit_rho <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned_rho$Q,
    phi    = tuned_rho$phi,
    solver = "cholesky"
  )
  fit_rho <- fastblm::constrain(fit_rho, C_aug_full)

  r_hat    <- fit_rho$posterior_mean[seq_len(p)]
  beta_hat <- fit_rho$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit_rho, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit_rho)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))

  results_sif_cv_blocked_rho_cv <- list(
    run_name              = "sif_cv_blocked_rho_cv",
    tags                  = list(tuning     = "cv_blocked",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 n_blocks   = max(blocked_folds)),
    timestamp             = Sys.time(),
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit_rho$sigma2e,
    phi                   = tuned_rho$phi,
    rho_opt               = rho_hat,
    cv_curve              = tuned_rho$history,
    fold_assignments      = blocked_folds,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )

  usethis::use_data(results_sif_cv_blocked_rho_cv, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned_rho$phi, fit_rho$sigma2e))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF blocked CV results saved to data/.")
