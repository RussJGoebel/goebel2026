# data-raw/run_sensitivity_sif_cv_blocked_weighted.R
#
# SIF sensitivity: weighted vs unweighted blocked CV with R_inv and free rho.
#
# Tests whether weighting CV scores by fold size (number of soundings per
# block) changes phi selection relative to equal-weight blocked CV.
# Both runs: water + RSR, R_inv, rho CV-tuned, spatially blocked folds.
#
# Requires the patched cv.R with weighted_folds parameter.
#
# Outputs:
#   results_sif_cv_blocked_weighted   -- blocked CV, weighted by fold size
#   results_sif_cv_blocked_unweighted -- blocked CV, equal-weight (baseline)

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

A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water
fine_grid_buffered <- d_shared$fine_grid_buffered

y_sif  <- goebel2026::soundings_augmented$SIF_757nm
R_inv  <- d_sif$R_inv

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

Q_fun_rho_aug <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  Q_a <- Matrix::forceSymmetric(
    Matrix::bdiag(Q, lambda_beta * Matrix::Diagonal(q))
  )
  list(Q = Q_a, log_det_Q = NULL)
}

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

blocked_folds <- d_sif$blocked_folds %||% goebel2026::setup_albedo$blocked_folds
if (is.null(blocked_folds)) stop("blocked_folds not found")

message(sprintf("Blocked folds: %d blocks", max(blocked_folds)))
message(sprintf("Block sizes: %s", paste(table(blocked_folds), collapse = " ")))

# ------------------------------------------------------------------------------
# 3. Shared fit helper
# ------------------------------------------------------------------------------

fit_cov <- function(tuned, rho_val) {
  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    R_inv  = R_inv,
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

  list(
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    cv_curve              = tuned$history,
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
# 5. Run 1: weighted blocked CV
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_cv_blocked_weighted")) {
  message("\n== SIF blocked CV: weighted by fold size, R_inv, rho tuned ==")

  tuned_w <- fastblm::tune_cv(
    y              = y_sif,
    A              = A_aug,
    Q_fun          = Q_fun_rho_aug,
    R_inv          = R_inv,
    theta_init     = c(rho = 0.9),
    lower          = c(rho = 0.5),
    upper          = c(rho = 0.999),
    k              = 10L,
    solver         = "cholesky",
    folds          = blocked_folds,
    constraint     = rsr_constraint,
    log_phi_lower  = log(0.001),
    weighted_folds = TRUE,
    parallel       = TRUE,
    verbose        = TRUE
  )

  rho_hat_w <- tuned_w$theta[["rho"]]
  message(sprintf("  weighted: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_w, tuned_w$phi, tuned_w$sigma2e))

  out_w <- fit_cov(tuned_w, rho_hat_w)

  results_sif_cv_blocked_weighted <- c(
    list(run_name  = "sif_cv_blocked_weighted",
         tags      = list(tuning     = "cv_blocked_weighted",
                          response   = "SIF_757nm",
                          covariates = "water",
                          constraint = "RSR",
                          W          = "queen",
                          rho        = rho_hat_w,
                          R_inv      = "sif_uncertainty",
                          n_blocks   = max(blocked_folds)),
         timestamp = Sys.time()),
    out_w
  )

  usethis::use_data(results_sif_cv_blocked_weighted, overwrite = TRUE)
  message(sprintf("  rho=%.4f  phi=%.4f", rho_hat_w, tuned_w$phi))
}

# ------------------------------------------------------------------------------
# 6. Run 2: unweighted blocked CV (equal-weight baseline)
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_cv_blocked_unweighted")) {
  message("\n== SIF blocked CV: equal-weight, R_inv, rho tuned ==")

  tuned_uw <- fastblm::tune_cv(
    y              = y_sif,
    A              = A_aug,
    Q_fun          = Q_fun_rho_aug,
    R_inv          = R_inv,
    theta_init     = c(rho = 0.9),
    lower          = c(rho = 0.5),
    upper          = c(rho = 0.999),
    k              = 10L,
    solver         = "cholesky",
    folds          = blocked_folds,
    constraint     = rsr_constraint,
    log_phi_lower  = log(0.001),
    weighted_folds = FALSE,
    parallel       = TRUE,
    verbose        = TRUE
  )

  rho_hat_uw <- tuned_uw$theta[["rho"]]
  message(sprintf("  unweighted: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_uw, tuned_uw$phi, tuned_uw$sigma2e))

  out_uw <- fit_cov(tuned_uw, rho_hat_uw)

  results_sif_cv_blocked_unweighted <- c(
    list(run_name  = "sif_cv_blocked_unweighted",
         tags      = list(tuning     = "cv_blocked_unweighted",
                          response   = "SIF_757nm",
                          covariates = "water",
                          constraint = "RSR",
                          W          = "queen",
                          rho        = rho_hat_uw,
                          R_inv      = "sif_uncertainty",
                          n_blocks   = max(blocked_folds)),
         timestamp = Sys.time()),
    out_uw
  )

  usethis::use_data(results_sif_cv_blocked_unweighted, overwrite = TRUE)
  message(sprintf("  rho=%.4f  phi=%.4f", rho_hat_uw, tuned_uw$phi))
}

# ------------------------------------------------------------------------------
# 7. Summary comparison
# ------------------------------------------------------------------------------

message("\n== Summary ==")
message(sprintf("  block sizes: %s", paste(table(blocked_folds), collapse = " ")))
if (exists("tuned_w"))
  message(sprintf("  weighted:   rho=%.4f  phi=%.4f  tau=%.3f",
                  tuned_w$theta[["rho"]], tuned_w$phi, 1/tuned_w$phi))
if (exists("tuned_uw"))
  message(sprintf("  unweighted: rho=%.4f  phi=%.4f  tau=%.3f",
                  tuned_uw$theta[["rho"]], tuned_uw$phi, 1/tuned_uw$phi))

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSIF blocked CV weighted vs unweighted results saved to data/.")
