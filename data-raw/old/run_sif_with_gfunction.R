# data-raw/run_sensitivity_sif_gA.R
#
# SIF sensitivity: ML vs CV vs blocked CV using g-weighted forward operator.
#
# Mirrors run_sensitivity_ml_sif.R and run_sensitivity_sif_cv_blocked.R but
# replaces A_flat with the Gaussian-weighted A_g (tau=1/3). Purpose: test
# whether the Gaussian forward model closes the gap between ML and CV,
# and whether blocked CV still collapses to tiny phi.
#
# All runs: queen adjacency, water + RSR (where applicable), rho=1 for CV,
# rho free for ML.
#
# Outputs:
#   results_sif_gA_cv          -- random CV, water + RSR, R_inv, rho=1
#   results_sif_gA_cv_norinv   -- random CV, water + RSR, no R_inv, rho=1
#   results_sif_gA_ml          -- ML, water, no RSR, rho free
#   results_sif_gA_ml_rho99    -- ML, water, no RSR, rho=0.99 fixed
#   results_sif_gA_cv_blocked  -- blocked CV, water + RSR, R_inv, rho=1

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
d_sif    <- goebel2026::setup_sif
d_gA     <- goebel2026::setup_g_A

# Swap A_flat for g-weighted A
A                  <- as(d_gA$A_g[["tau_0.333"]], "dgCMatrix")
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water
fine_grid_buffered <- d_shared$fine_grid_buffered

y_sif  <- goebel2026::soundings_augmented$SIF_757nm
R_inv  <- d_sif$R_inv

blocked_folds <- d_sif$blocked_folds %||%
  goebel2026::setup_albedo$blocked_folds

`%||%` <- function(x, y) if (!is.null(x)) x else y

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

message(sprintf("A_g dim: %d x %d  (target pixels: %d)",
                nrow(A), ncol(A), length(target_idx)))

# ------------------------------------------------------------------------------
# 2. Q functions
# ------------------------------------------------------------------------------

# Intrinsic SAR (rho=1) -- for CV runs
Q_fun_fixed_aug <- local({
  S   <- Matrix::Diagonal(nrow(W_queen)) - W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  Q_a <- Matrix::forceSymmetric(
    Matrix::bdiag(Q, lambda_beta * Matrix::Diagonal(q))
  )
  function(theta) list(Q = Q_a, log_det_Q = NULL)
})

# Rho-tunable -- for ML runs (spatial block only; tune_ml augments internally)
Q_fun_rho <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  list(Q = Q, log_det_Q = NULL)
}

# Fixed rho=0.99 -- for ML fixed-rho run
Q_fun_fixed_99 <- local({
  S <- Matrix::Diagonal(nrow(W_queen)) - 0.99 * W_queen
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
})

# RSR constraint per fold
rsr_constraint <- local({
  force(X_obs_water); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs_water[train_idx, , drop = FALSE]) %*% A_train)
    C_zeros   <- matrix(0, nrow = q, ncol = q)
    cbind(C_spatial, C_zeros)
  }
})

C_aug_full <- cbind(
  as.matrix(t(X_obs_water) %*% A),
  matrix(0, nrow = q, ncol = q)
)

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_sif), 10L)

# ------------------------------------------------------------------------------
# 3. Fit helpers
# ------------------------------------------------------------------------------

fit_cov_rsr <- function(tuned, rho_val, R_inv_fit = NULL) {
  Q_aug <- Matrix::forceSymmetric(
    Matrix::bdiag(tuned$Q[seq_len(p), seq_len(p)],
                  lambda_beta * Matrix::Diagonal(q))
  )
  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = tuned$Q, phi = tuned$phi,
    R_inv = R_inv_fit, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

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

fit_cov_ml <- function(tuned, rho_val) {
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

# --- 1. Random CV, water + RSR, R_inv, rho=1 ----------------------------------

if (!.should_skip("results_sif_gA_cv")) {
  message("\n== SIF g-A CV: water + RSR, R_inv, rho=1 ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    R_inv      = R_inv,
    theta_init = numeric(0),
    k          = 10L,
    solver     = "cholesky",
    folds      = fold_assignments,
    constraint = rsr_constraint,
    parallel   = TRUE,
    verbose    = TRUE
  )
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))
  out <- fit_cov_rsr(tuned, rho_val = 1, R_inv_fit = R_inv)
  results_sif_gA_cv <- c(
    list(run_name = "sif_gA_cv",
         tags     = list(tuning="cv", response="SIF_757nm", covariates="water",
                         constraint="RSR", W="queen", rho=1, A="g_tau0.333",
                         R_inv="sif_uncertainty"),
         timestamp = Sys.time(), cv_history = tuned$history),
    out)
  usethis::use_data(results_sif_gA_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f", tuned$phi))
}

# --- 2. Random CV, water + RSR, no R_inv, rho=1 -------------------------------

if (!.should_skip("results_sif_gA_cv_norinv")) {
  message("\n== SIF g-A CV: water + RSR, no R_inv, rho=1 ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    theta_init = numeric(0),
    k          = 10L,
    solver     = "cholesky",
    folds      = fold_assignments,
    constraint = rsr_constraint,
    parallel   = TRUE,
    verbose    = TRUE
  )
  message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))
  out <- fit_cov_rsr(tuned, rho_val = 1, R_inv_fit = NULL)
  results_sif_gA_cv_norinv <- c(
    list(run_name = "sif_gA_cv_norinv",
         tags     = list(tuning="cv", response="SIF_757nm", covariates="water",
                         constraint="RSR", W="queen", rho=1, A="g_tau0.333",
                         R_inv=NULL),
         timestamp = Sys.time(), cv_history = tuned$history),
    out)
  usethis::use_data(results_sif_gA_cv_norinv, overwrite = TRUE)
  message(sprintf("  phi=%.4f", tuned$phi))
}

# --- 3. ML, water, no RSR, rho free -------------------------------------------

if (!.should_skip("results_sif_gA_ml")) {
  message("\n== SIF g-A ML: water, rho free ==")

  tuned <- fastblm::tune_ml(
    y          = y_sif,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = X_obs_water,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )
  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))
  out <- fit_cov_ml(tuned, rho_hat)
  results_sif_gA_ml <- c(
    list(run_name  = "sif_gA_ml",
         tags      = list(tuning="ml", response="SIF_757nm", covariates="water",
                          constraint="none", W="queen", rho=rho_hat,
                          A="g_tau0.333"),
         timestamp  = Sys.time(), ml_history = tuned$history),
    out)
  usethis::use_data(results_sif_gA_ml, overwrite = TRUE)
  message(sprintf("  rho=%.4f  phi=%.4f", rho_hat, tuned$phi))
}

# --- 4. ML, water, no RSR, rho=0.99 fixed -------------------------------------

if (!.should_skip("results_sif_gA_ml_rho99")) {
  message("\n== SIF g-A ML: water, rho=0.99 fixed ==")

  tuned <- fastblm::tune_ml(
    y       = y_sif,
    A       = A,
    Q_fun   = Q_fun_fixed_99,
    X_fixed = X_obs_water,
    verbose = TRUE
  )
  message(sprintf("  rho=0.99 (fixed)  phi=%.4f  sigma2e=%.4g",
                  tuned$phi, tuned$sigma2e))
  out <- fit_cov_ml(tuned, rho_val = 0.99)
  results_sif_gA_ml_rho99 <- c(
    list(run_name  = "sif_gA_ml_rho99",
         tags      = list(tuning="ml", response="SIF_757nm", covariates="water",
                          constraint="none", W="queen", rho=0.99,
                          A="g_tau0.333"),
         timestamp  = Sys.time(), ml_history = tuned$history),
    out)
  usethis::use_data(results_sif_gA_ml_rho99, overwrite = TRUE)
  message(sprintf("  phi=%.4f", tuned$phi))
}

# --- 5. Blocked CV, water + RSR, R_inv, rho=1 ---------------------------------

if (!.should_skip("results_sif_gA_cv_blocked")) {
  message("\n== SIF g-A blocked CV: water + RSR, R_inv, rho=1 ==")

  if (is.null(blocked_folds)) {
    message("  blocked_folds not found -- skipping.")
  } else {
    message(sprintf("  %d blocks, sizes: %s",
                    max(blocked_folds),
                    paste(table(blocked_folds), collapse = " ")))

    tuned <- fastblm::tune_cv(
      y             = y_sif,
      A             = A_aug,
      Q_fun         = Q_fun_fixed_aug,
      R_inv         = R_inv,
      theta_init    = numeric(0),
      k             = 10L,
      solver        = "cholesky",
      folds         = blocked_folds,
      constraint    = rsr_constraint,
      log_phi_lower = log(0.001),
      parallel      = TRUE,
      verbose       = TRUE
    )
    message(sprintf("  phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))
    out <- fit_cov_rsr(tuned, rho_val = 1, R_inv_fit = R_inv)
    results_sif_gA_cv_blocked <- c(
      list(run_name = "sif_gA_cv_blocked",
           tags     = list(tuning="cv_blocked", response="SIF_757nm",
                           covariates="water", constraint="RSR", W="queen",
                           rho=1, A="g_tau0.333", R_inv="sif_uncertainty",
                           n_blocks=max(blocked_folds)),
           timestamp = Sys.time(), cv_history = tuned$history),
      out)
    usethis::use_data(results_sif_gA_cv_blocked, overwrite = TRUE)
    message(sprintf("  phi=%.4f", tuned$phi))
  }
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF g-A sensitivity results saved to data/.")
