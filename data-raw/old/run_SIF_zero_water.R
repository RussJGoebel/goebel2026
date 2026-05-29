# data-raw/run_sif_water_zero.R
#
# SIF downscaling run with heterogeneous SAR prior (water-zero).
#
# Collaborators' prior belief: water pixels have zero SIF.
# Implemented via a heterogeneous SAR prior where the innovation variance
# at each pixel is proportional to its land fraction:
#
#   Q_rho = (I - rho*W)' * diag(1/land_frac) * (I - rho*W)
#
# No water covariate -- plain A (not A_aug).
# rho tuned via CV golden section, phi profiled at each rho.
# No ground truth -- real data run, posterior mean and SE only.
#
# Outputs saved to data/ via usethis::use_data():
#   results_sif_water_zero

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)
message(sprintf("Using %d workers", parallel::detectCores() - 1L))

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_latent_water     <- d_shared$X_latent_water

y_sif <- goebel2026::soundings_augmented$SIF_757nm

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

land_frac      <- 1 - as.numeric(X_latent_water[, 2])
land_frac_safe <- pmax(land_frac, 1e-3)

message(sprintf("Land fraction summary: min=%.4f  mean=%.4f  max=%.4f",
                min(land_frac), mean(land_frac), max(land_frac)))
message(sprintf("Pure water pixels (land_frac < 0.01): %d of %d",
                sum(land_frac < 0.01), p))

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_sif), 10L)

# ------------------------------------------------------------------------------
# 3. Q_fun -- heterogeneous SAR, explicit sparse matrix for Cholesky solver
# ------------------------------------------------------------------------------

Q_fun_sif_water_zero <- function(theta) {
  rho <- min(max(theta[["rho"]], 0.01), 0.999)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  D_inv <- Matrix::Diagonal(x = land_frac_safe)
  Q_rho <- Matrix::drop0(
    Matrix::forceSymmetric(Matrix::crossprod(S_rho, D_inv) %*% S_rho)
  )

  list(Q = Q_rho, Q_rho = Q_rho, log_det_Q = NULL, precond = NULL)
}

# ------------------------------------------------------------------------------
# 4. Sanity check timing
# ------------------------------------------------------------------------------

message("\nTiming sanity check (sequential, k=10, phi=100)...")
prior_check <- Q_fun_sif_water_zero(c(rho = 0.9))
t_check <- system.time({
  cv_check <- fastblm:::.eval_cv(
    y_sif, A, Q_fun_sif_water_zero, c(rho = 0.9), prior_check,
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

if (!.should_skip("results_sif_water_zero")) {
  message("\n== SIF water-zero: heterogeneous SAR, no water covariate ==")

  tuned <- fastblm::tune_cv(
    y           = y_sif,
    A           = A,
    Q_fun       = Q_fun_sif_water_zero,
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
  prior_final <- Q_fun_sif_water_zero(c(rho = rho_hat))

  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A,
    Q      = prior_final$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )

  mu <- as.numeric(fit$posterior_mean)

  A_pred <- Matrix::Diagonal(p)
  se     <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)

  mu_t <- mu[target_idx]
  se_t <- se[target_idx]
  ns_t <- n_soundings_per_pixel[target_idx]

  results_sif_water_zero <- list(
    run_name              = "sif_water_zero",
    tags                  = list(tuning     = "cv",
                                 response   = "SIF_757nm",
                                 covariates = "none",
                                 constraint = "heterogeneous_SAR",
                                 W          = "queen",
                                 rho        = rho_hat),
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    land_frac             = land_frac[target_idx],
    cv_curve              = tuned$history,
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_sif_water_zero, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, fit$sigma2e))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF water-zero results saved to data/.")
