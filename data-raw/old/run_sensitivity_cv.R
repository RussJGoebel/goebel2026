# data-raw/run_sensitivity_cv.R
#
# Supplement S5: Spatially blocked cross-validation sensitivity.
# Compares random CV folds (baseline) vs spatially blocked folds (k-means
# clusters of sounding centroids). Both use: water + RSR + rho=1.
# Baseline result (random CV) is results_water_rsr_rho1 from run_main_results.R.
#
# Expected finding: blocked CV selects similar tau but may undersmooth due
# to target-mode observation density concentrating overlap in the center --
# blocking mostly withholds the overlap-rich region, causing the CV objective
# to favour less smoothing.
#
# Outputs saved to data/ via usethis::use_data():
#   results_cv_blocked   -- spatially blocked folds

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- FALSE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_alb          <- d_albedo$y
y_latent_true  <- d_albedo$y_latent_true
blocked_folds  <- d_albedo$blocked_folds

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

message(sprintf("Blocked folds: %d blocks, sizes: %s",
                max(blocked_folds),
                paste(table(blocked_folds), collapse = " ")))

# ------------------------------------------------------------------------------
# 2. Shared model components
# ------------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_Q_fun_fixed <- function(W) {
  S <- Matrix::Diagonal(nrow(W)) - W
  Q <- Matrix::forceSymmetric(Matrix::crossprod(S))
  Q <- Matrix::drop0(Q)
  function(theta) list(Q = Q, log_det_Q = NULL)
}

make_Q_fun_aug <- function(Q_fun_spatial, q, lambda) {
  force(Q_fun_spatial); force(q); force(lambda)
  function(theta) {
    sp  <- Q_fun_spatial(theta)
    Q_a <- Matrix::forceSymmetric(
      Matrix::bdiag(sp$Q, lambda * Matrix::Diagonal(q))
    )
    list(Q = Q_a, log_det_Q = NULL)
  }
}

make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    C_zeros   <- matrix(0, nrow = q, ncol = q)
    cbind(C_spatial, C_zeros)
  }
}

Q_fun_fixed     <- make_Q_fun_fixed(W_queen)
Q_fun_fixed_aug <- make_Q_fun_aug(Q_fun_fixed, q, lambda_beta)
A_aug           <- as(cbind(A, X_obs_water), "dgCMatrix")
rsr_constraint  <- make_rsr_constraint(X_obs_water, p, q)
C_spatial_full  <- as.matrix(t(X_obs_water) %*% A)
C_zeros_full    <- matrix(0, nrow = q, ncol = q)
C_aug_full      <- cbind(C_spatial_full, C_zeros_full)

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

compute_outputs <- function(result, run_name, tags) {
  fit   <- result$fit
  tuned <- result$tuned

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred)
  se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t
  ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }

  resid  <- mu_t - tr_t
  rmse   <- sqrt(mean(resid^2, na.rm = TRUE))
  ss_res <- sum(resid^2, na.rm = TRUE)
  ss_tot <- sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)
  r2     <- 1 - ss_res / ss_tot

  list(
    run_name              = run_name,
    tags                  = tags,
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = 1,
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

# ------------------------------------------------------------------------------
# 3. Blocked CV run
# ------------------------------------------------------------------------------

if (!.should_skip("results_cv_blocked")) {
  message("\n== S5: Spatially blocked CV ==")
  message(sprintf("  %d blocks, sizes: %s",
                  max(blocked_folds),
                  paste(table(blocked_folds), collapse = " ")))

  tuned_blocked <- fastblm::tune_cv(
    y          = y_alb,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    theta_init = numeric(0),
    folds      = blocked_folds,
    constraint = rsr_constraint,
    parallel   = TRUE,
    verbose    = TRUE
  )

  fit_blocked <- fastblm::fit_fastblm(
    y      = y_alb,
    A      = A_aug,
    Q      = tuned_blocked$Q,
    phi    = tuned_blocked$phi,
    solver = "cholesky"
  )
  fit_blocked <- fastblm::constrain(fit_blocked, C_aug_full)

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g",
                  tuned_blocked$phi, tuned_blocked$sigma2e))

  out_blocked <- compute_outputs(
    list(fit = fit_blocked, tuned = tuned_blocked),
    "cv_blocked",
    tags = list(tuning = "cv_blocked", covariates = "water",
                constraint = "RSR", W = "queen", rho = 1,
                n_blocks = max(blocked_folds))
  )
  results_cv_blocked <- out_blocked
  usethis::use_data(results_cv_blocked, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  phi=%.4f",
                  out_blocked$rmse, out_blocked$r2,
                  out_blocked$coverage_95_obs, out_blocked$phi))
  message(sprintf("  phi comparison -- random CV baseline: %.4f  blocked: %.4f",
                  goebel2026::results_water_rsr_rho1$phi, tuned_blocked$phi))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nS5 blocked CV sensitivity results saved to data/.")
