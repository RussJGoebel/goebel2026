# data-raw/run_sif_rsr_rinv.R
#
# SIF downscaling run with hard RSR, rho CV-tuned, heteroskedastic R_inv.
#
# Same as run_sif_rsr.R but with R_inv = diag(1/sigma_eps^2) using
# SIF_Uncertainty_757nm as per-sounding standard errors.
#
# Outputs saved to data/ via usethis::use_data():
#   results_sif_rsr_rinv

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multicore, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_sif     <- goebel2026::soundings_augmented$SIF_757nm
sigma_eps <- goebel2026::soundings_augmented$SIF_Uncertainty_757nm
R_inv     <- Matrix::Diagonal(x = 1 / sigma_eps^2)

message(sprintf("Noise weights: min_sigma=%.4f  median_sigma=%.4f  max_sigma=%.4f",
                min(sigma_eps), median(sigma_eps), max(sigma_eps)))

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_sif), 10L)

# ------------------------------------------------------------------------------
# 3. Q_fun -- plain SAR + beta ridge, 1-arg
# ------------------------------------------------------------------------------

Q_fun_sif <- function(theta) {
  rho <- min(max(theta[["rho"]], 0.01), 0.999)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q_rho <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))

  # Augmented precision: SAR block + beta ridge
  Q_aug <- Matrix::drop0(Matrix::forceSymmetric(
    Matrix::bdiag(Q_rho, lambda_beta * Matrix::Diagonal(q))
  ))

  list(Q = Q_aug, Q_rho = Q_rho, log_det_Q = NULL, precond = NULL)
}

# ------------------------------------------------------------------------------
# 4. Fold-specific RSR constraint
#
# C is q x (p+q): spatial block [t(X_obs_train) %*% A_train] plus zero
# padding for the beta columns. The constraint only applies to the spatial
# field r, not to beta.
# ------------------------------------------------------------------------------

constraint_rsr <- function(train_idx, A_train) {
  A_sp_train <- A_train[, seq_len(p), drop = FALSE]
  X_tr       <- A_train[, p + seq_len(q), drop = FALSE]
  C_spatial  <- as.matrix(Matrix::crossprod(X_tr, A_sp_train))  # q x p
  cbind(C_spatial, matrix(0, nrow = q, ncol = q))               # q x (p+q)
}

# Full-data constraint for final fit
C_full <- {
  C_spatial_full <- as.matrix(Matrix::crossprod(X_obs_water, A))  # q x p
  cbind(C_spatial_full, matrix(0, nrow = q, ncol = q))            # q x (p+q)
}

# ------------------------------------------------------------------------------
# 5. Sanity check timing
# ------------------------------------------------------------------------------

message("\nTiming sanity check (sequential, k=10, phi=100)...")
prior_check <- Q_fun_sif(c(rho = 0.9))
fold_C_list <- fastblm:::.precompute_fold_constraints(
  fastblm:::.make_constraint_fn(constraint_rsr, p + q), A_aug, fold_assignments
)
t_check <- system.time({
  cv_check <- fastblm:::.eval_cv(
    y_sif, A_aug, Q_fun_sif, c(rho = 0.9), prior_check,
    100, fold_assignments,
    fastblm:::.make_score_fn("mse"),
    "cholesky", 1e-6, NULL, NULL,
    fold_C_list = fold_C_list,
    precond_fun = NULL,
    parallel    = FALSE
  )
})
message(sprintf("  %.2fs total  %.2fs/fold  cv_mse=%.4e",
                t_check["elapsed"], t_check["elapsed"] / 10, cv_check))

# ------------------------------------------------------------------------------
# 6. Run
# ------------------------------------------------------------------------------

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name))
    return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error   = function(e) NULL, warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}

if (!.should_skip("results_sif_rsr_rinv")) {
  message("\n== SIF RSR: rho CV-tuned, hard RSR, R_inv ==")

  tuned <- fastblm::tune_cv(
    y           = y_sif,
    A           = A_aug,
    Q_fun       = Q_fun_sif,
    R_inv       = R_inv,
    theta_init  = c(rho = 0.9),
    lower       = c(rho = 0.5),
    upper       = c(rho = 0.999),
    k           = 10L,
    solver      = "cholesky",
    constraint  = constraint_rsr,
    precond_fun = NULL,
    parallel    = TRUE,
    verbose     = TRUE,
    folds       = fold_assignments
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  # Final fit with full-data constraint
  prior_final <- Q_fun_sif(c(rho = rho_hat))

  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = prior_final$Q,
    phi    = tuned$phi,
    R_inv  = R_inv,
    solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_full)

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)

  A_pred  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se      <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]

  results_sif_rsr_rinv <- list(
    run_name              = "sif_rsr_rinv",
    tags                  = list(tuning     = "cv",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "hard_RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 R_inv      = "SIF_Uncertainty_757nm"),
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    cv_curve              = tuned$history,
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_sif_rsr_rinv, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g  beta_water=%.4f",
                  rho_hat, tuned$phi, fit$sigma2e, beta_hat[2]))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF RSR R_inv results saved to data/.")
