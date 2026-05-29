# data-raw/run_water_zero.R
#
# "Water zero" sensitivity run.
#
# Collaborators' prior belief: water pixels have zero SIF.
# Implemented via a heterogeneous SAR prior where the innovation variance
# at each pixel is proportional to its land fraction:
#
#   Q_rho = (I - rho*W)' * diag(1/land_frac) * (I - rho*W)
#
# Pure land pixels (land_frac=1): full SAR variance.
# Mixed pixels (land_frac in (0,1)): variance proportional to land_frac.
# Pure water pixels (land_frac=0): clamped to land_frac=1e-3, effectively
#   pinned to zero.
#
# No water covariate -- plain A (not A_aug).
# rho tuned via CV golden section (1D optimize, no gradient issues).
# phi profiled out at each rho.
#
# Outputs saved to data/ via usethis::use_data():
#   results_water_zero  -- heterogeneous SAR, no water covariate

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
d_albedo <- goebel2026::setup_albedo

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

# Land fraction: 1 = pure land, 0 = pure water
# Clamped at 1e-3 to avoid infinite precision at pure water pixels
land_frac      <- 1 - as.numeric(X_latent_water[, 2])
land_frac_safe <- pmax(land_frac, 1e-3)

message(sprintf("Land fraction summary: min=%.4f  mean=%.4f  max=%.4f",
                min(land_frac), mean(land_frac), max(land_frac)))
message(sprintf("Pure water pixels (land_frac < 0.01): %d of %d",
                sum(land_frac < 0.01), p))

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_alb), 10L)

# ------------------------------------------------------------------------------
# 3. Q_fun -- heterogeneous SAR, explicit sparse matrix for Cholesky solver
#
# Q = (I - rho*W)' * diag(1/land_frac) * (I - rho*W)
#
# Returned as an explicit sparse dgCMatrix so tune_cv can use the Cholesky
# path (p x p sparse Cholesky), which is exact and fast when A'A is sparse.
# ------------------------------------------------------------------------------

Q_fun_water_zero <- function(theta) {
  rho <- theta[["rho"]]
  rho <- min(max(rho, 0.01), 0.999)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  D_inv <- Matrix::Diagonal(x = 1 / land_frac_safe)
  Q_rho <- Matrix::drop0(
    Matrix::forceSymmetric(Matrix::crossprod(S_rho, D_inv) %*% S_rho)
  )

  list(Q = Q_rho, Q_rho = Q_rho, log_det_Q = NULL, precond = NULL)
}

# ------------------------------------------------------------------------------
# 4. Sanity check
# ------------------------------------------------------------------------------

message("\nSanity check (rho=0.9, phi=100, sequential)...")
prior_check <- Q_fun_water_zero(c(rho = 0.9))
t_check <- system.time({
  cv_check <- fastblm:::.eval_cv(
    y_alb, A, Q_fun_water_zero, c(rho = 0.9), prior_check,
    100, fold_assignments,
    fastblm:::.make_score_fn("mse"),
    "cholesky", 1e-6, NULL, NULL,
    fold_C_list = NULL, precond_fun = NULL, parallel = FALSE
  )
})
message(sprintf("  %.2fs total  %.2fs/fold  cv_mse=%.4e",
                t_check["elapsed"], t_check["elapsed"] / 10, cv_check))

# ------------------------------------------------------------------------------
# 5. Run
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

if (!.should_skip("results_water_zero")) {
  message("\n== Water-zero run: heterogeneous SAR, no water covariate ==")

  tuned <- fastblm::tune_cv(
    y           = y_alb,
    A           = A,
    Q_fun       = Q_fun_water_zero,
    theta_init  = c(rho = 0.9),
    lower       = c(rho = 0.5),
    upper       = c(rho = 0.999),
    k           = 10L,
    solver      = "cholesky",
    constraint  = NULL,
    precond_fun = NULL,
    parallel    = TRUE,
    verbose     = TRUE,
    folds       = fold_assignments
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  # Final fit
  prior_final <- Q_fun_water_zero(c(rho = rho_hat))

  fit <- fastblm::fit_fastblm(
    y      = y_alb,
    A      = A,
    Q      = prior_final$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )

  mu <- as.numeric(fit$posterior_mean)

  A_pred <- Matrix::Diagonal(p)
  se     <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)

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

  results_water_zero <- list(
    run_name              = "water_zero",
    tags                  = list(tuning     = "cv",
                                 covariates = "none",
                                 constraint = "heterogeneous_SAR",
                                 W          = "queen",
                                 rho        = rho_hat),
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    land_frac             = land_frac[target_idx],
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_water_zero, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat, rmse, r2, results_water_zero$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nWater-zero results saved to data/.")
