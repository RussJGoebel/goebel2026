# data-raw/run_150m_albedo.R
#
# 150m albedo validation experiment.
# Compares:
#   1. Bayesian downscaling with correct A (uniform intersection weights)
#   2. Bayesian downscaling with centroid A (point-support approximation)
#   3. Universal kriging (point support at sounding centroids)
#
# Uses PCG throughout (Cholesky infeasible at 150m due to ~80GB fill-in).
# All three methods use water proportion covariate.
# Methods 1+2 use RSR constraint and rho CV-tuned.
#
# Outputs saved via usethis::use_data():
#   results_water_rsr_rho_cv_150m     -- correct A, Bayesian
#   results_albedo_gA_centroid_150m   -- centroid A, Bayesian
#   results_kriging_albedo_150m       -- universal kriging

library(fastblm)
library(goebel2026)
library(spatintegrate)
library(Matrix)
library(sf)
library(gstat)
library(sp)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multisession())

# ------------------------------------------------------------------------------
# 1. Load setup
# ------------------------------------------------------------------------------

message("== Loading 150m setup ==")

d_shared <- goebel2026::setup_shared_150m
d_albedo <- goebel2026::setup_albedo_150m

fine_grid_buffered <- d_shared$fine_grid_buffered
A_flat             <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A_flat)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A_flat > 0))

message(sprintf("  Grid cells total: %d  target: %d  soundings: %d",
                nrow(fine_grid_buffered), length(target_idx), nrow(A_flat)))

# ------------------------------------------------------------------------------
# 2. Shared model components
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A_flat, X_obs_water), "dgCMatrix")

# Q_fun: fixed rho=1 (intrinsic SAR), augmented with beta prior
make_Q_fun_rho1_aug <- function(W, q, lambda) {
  force(W); force(q); force(lambda)
  S    <- Matrix::Diagonal(nrow(W)) - W
  Q_sp <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  Q_a  <- Matrix::drop0(Matrix::forceSymmetric(
    Matrix::bdiag(Q_sp, lambda * Matrix::Diagonal(q))
  ))
  Q_fixed <- Q_a
  function(theta) list(Q = Q_fixed, log_det_Q = NULL)
}

Q_fun_rho1_aug <- make_Q_fun_rho1_aug(W_queen, q, lambda_beta)

# Random CV folds
set.seed(2026L)
folds_random <- sample(rep_len(1:10, length(y_alb)))

# RSR constraint: per-fold version
make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    cbind(C_spatial, matrix(0, nrow = q, ncol = q))
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

# Full-data RSR constraint for final fit
C_aug_full <- cbind(
  as.matrix(t(X_obs_water) %*% A_flat),
  matrix(0, nrow = q, ncol = q)
)

# ------------------------------------------------------------------------------
# Helper: build augmented Q from tuned spatial Q
# tune_cv returns only the spatial Q block; we need bdiag(Q_spatial, lambda*I_q)
# for fit_fastblm when A is the augmented [A_flat | X_obs_water] matrix.
# Without this, apply_K operates on a (p+q)-vector with a p x p Q, producing
# a malformed PCG system that never converges correctly.
# ------------------------------------------------------------------------------
make_Q_aug <- function(Q_spatial, q, lambda) {
  Matrix::drop0(Matrix::forceSymmetric(
    Matrix::bdiag(Q_spatial, lambda * Matrix::Diagonal(q))
  ))
}

# ------------------------------------------------------------------------------
# Helper: sparse Cholesky preconditioner for K = A'A + (1/phi)*Q_aug
#
# Uses Cholesky of (Q_aug + eps*I) as preconditioner. This is highly effective
# because:
#   - For unobserved pixels, K ≈ (1/phi)*Q so Q^{-1} is nearly exact there
#   - The diagonal fallback handles observed pixels adequately
#   - eps regularisation is needed because the intrinsic SAR (rho=1) makes
#     the spatial block of Q singular (rank p-1, null space = constant vector)
#
# The Cholesky is computed once and reused across all PCG iterations.
# Cost per PCG step: one sparse triangular solve -- O(nnz(Q)), very cheap.
# ------------------------------------------------------------------------------
make_chol_precond <- function(Q_aug, eps = 1e-6) {
  Q_reg  <- Q_aug + eps * Matrix::Diagonal(nrow(Q_aug))
  CQ     <- Matrix::Cholesky(Q_reg, LDL = FALSE, perm = TRUE)
  function(v) as.numeric(Matrix::solve(CQ, v))
}

# ------------------------------------------------------------------------------
# 3. Output helper
# ------------------------------------------------------------------------------

compute_outputs_albedo <- function(fit, tuned, run_name, tags) {
  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")

  message("  Computing posterior SE (PCG+Lanczos, n_probes=200)...")
  se      <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

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

  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  beta=[%.4f, %.4f]",
                  rmse, r2, coverage(which(ns_t >= 1L)),
                  beta_hat[1], beta_hat[2]))

  list(
    run_name              = run_name,
    tags                  = tags,
    timestamp             = Sys.time(),
    resolution_m          = 150L,
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = 1.0,
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
# 4. Method 1: Correct A (uniform intersection weights)
# ------------------------------------------------------------------------------

if (FORCE_RERUN || !exists("results_water_rsr_rho_cv_150m")) {
  message("\n== 1. Bayesian (correct A, 150m) ==")

  tuned <- fastblm::tune_cv(
    y          = y_alb,
    A          = A_aug,
    Q_fun      = Q_fun_rho1_aug,
    theta_init = numeric(0),
    k          = 10L,
    constraint = rsr_constraint,
    folds      = folds_random,
    solver     = "pcg",
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat <- 1.0
  message(sprintf("  tuned: rho=1 (fixed)  phi=%.4f", tuned$phi))

  # FIX: tune_cv returns only the spatial Q block (p x p). The design matrix
  # A_aug is (p+q)-wide, so fit_fastblm needs the full augmented Q.
  Q_aug_final <- make_Q_aug(tuned$Q, q, lambda_beta)

  # Preconditioner: sparse Cholesky of Q_aug (with small ridge for rho=1 singularity).
  # Factored once here, reused across all PCG iterations inside fit_fastblm.
  message("  Building Cholesky preconditioner...")
  precond_1 <- make_chol_precond(Q_aug_final)

  fit <- fastblm::fit_fastblm(
    y           = y_alb,
    A           = A_aug,
    Q           = Q_aug_final,   # was: tuned$Q  (p x p only -- dimension mismatch)
    phi         = tuned$phi,
    solver      = "pcg",
    pcg_precond = precond_1
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_water_rsr_rho_cv_150m <- compute_outputs_albedo(
    fit, tuned, "water_rsr_rho_cv_150m",
    tags = list(resolution = 150L, tuning = "cv", covariates = "water",
                constraint = "RSR", W = "queen", A = "uniform",
                solver = "pcg")
  )
  usethis::use_data(results_water_rsr_rho_cv_150m, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 5. Method 2: Centroid A (point-support approximation)
# ------------------------------------------------------------------------------

if (FORCE_RERUN || !exists("results_albedo_gA_centroid_150m")) {
  message("\n== 2. Bayesian (centroid A, 150m) ==")

  # Build centroid A: each sounding contributes only to the cell
  # containing its centroid, with weight 1
  message("  Building centroid A...")
  snd_cents  <- sf::st_centroid(d_shared$soundings_proj)
  grid_geom  <- sf::st_geometry(fine_grid_buffered)
  hits       <- sf::st_within(snd_cents, fine_grid_buffered, sparse = TRUE)

  i_idx <- which(lengths(hits) > 0)
  j_idx <- unlist(hits[i_idx])
  A_centroid <- Matrix::sparseMatrix(
    i    = i_idx,
    j    = j_idx,
    x    = 1.0,
    dims = c(nrow(A_flat), ncol(A_flat))
  )
  A_centroid <- as(A_centroid, "dgCMatrix")
  message(sprintf("  A_centroid: %d x %d  nnz=%d  (%.1f%% soundings matched)",
                  nrow(A_centroid), ncol(A_centroid),
                  Matrix::nnzero(A_centroid),
                  100 * length(i_idx) / nrow(A_flat)))

  A_aug_cent <- as(cbind(A_centroid, X_obs_water), "dgCMatrix")

  C_aug_cent <- cbind(
    as.matrix(t(X_obs_water) %*% A_centroid),
    matrix(0, nrow = q, ncol = q)
  )

  rsr_constraint_cent <- make_rsr_constraint(X_obs_water, p, q)

  Q_fun_cent <- make_Q_fun_rho1_aug(W_queen, q, lambda_beta)

  tuned_cent <- fastblm::tune_cv(
    y          = y_alb,
    A          = A_aug_cent,
    Q_fun      = Q_fun_cent,
    theta_init = numeric(0),
    k          = 10L,
    constraint = rsr_constraint_cent,
    folds      = folds_random,
    solver     = "pcg",
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat_cent <- 1.0
  message(sprintf("  tuned: rho=1 (fixed)  phi=%.4f", tuned_cent$phi))

  # FIX: same issue -- tuned_cent$Q is spatial-only, need full augmented Q
  Q_aug_cent_final <- make_Q_aug(tuned_cent$Q, q, lambda_beta)

  message("  Building Cholesky preconditioner...")
  precond_2 <- make_chol_precond(Q_aug_cent_final)

  fit_cent <- fastblm::fit_fastblm(
    y           = y_alb,
    A           = A_aug_cent,
    Q           = Q_aug_cent_final,  # was: tuned_cent$Q  (p x p only)
    phi         = tuned_cent$phi,
    solver      = "pcg",
    pcg_precond = precond_2
  )
  fit_cent <- fastblm::constrain(fit_cent, C_aug_cent)

  # Use same compute_outputs but with centroid A
  r_hat    <- fit_cent$posterior_mean[seq_len(p)]
  beta_hat <- fit_cent$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")

  message("  Computing posterior SE...")
  se      <- fastblm::posterior_se(fit_cent, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit_cent, n_probes = 200L)[p + seq_len(q)]

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
  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rmse, r2, coverage(which(ns_t >= 1L))))

  results_albedo_gA_centroid_150m <- list(
    run_name              = "albedo_gA_centroid_150m",
    tags                  = list(resolution = 150L, tuning = "cv",
                                 covariates = "water", constraint = "RSR",
                                 W = "queen", A = "centroid", solver = "pcg"),
    timestamp             = Sys.time(),
    resolution_m          = 150L,
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit_cent$sigma2e,
    phi                   = tuned_cent$phi,
    rho_opt               = 1.0,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_albedo_gA_centroid_150m, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 6. Method 3: Universal kriging at 150m prediction locations
# ------------------------------------------------------------------------------

if (FORCE_RERUN || !exists("results_kriging_albedo_150m")) {
  message("\n== 3. Universal kriging (150m) ==")

  sounding_centroids <- sf::st_centroid(d_shared$soundings_proj)
  water_obs          <- X_obs_water[, 2]
  target_sf          <- fine_grid_buffered[target_idx, ]
  target_cents       <- sf::st_centroid(target_sf)
  water_target       <- target_sf$proportion_water

  sounding_sp <- as(sounding_centroids, "Spatial")
  target_sp   <- as(target_cents,       "Spatial")

  obs_df <- data.frame(y = y_alb, water = water_obs)
  coordinates(obs_df) <- coordinates(sounding_sp)
  proj4string(obs_df)  <- proj4string(sounding_sp)

  pred_df <- data.frame(water = water_target)
  coordinates(pred_df) <- coordinates(target_sp)
  proj4string(pred_df)  <- proj4string(target_sp)

  message("  Fitting empirical variogram...")
  vgm_emp <- gstat::variogram(y ~ water, data = obs_df,
                              cutoff = 15000, width = 500)

  psill_init  <- var(y_alb, na.rm = TRUE) * 0.8
  range_init  <- 3000
  nugget_init <- var(y_alb, na.rm = TRUE) * 0.1

  vgm_fit <- gstat::fit.variogram(
    vgm_emp,
    model = gstat::vgm(psill  = psill_init,
                       model  = "Exp",
                       range  = range_init,
                       nugget = nugget_init)
  )
  message(sprintf("  Variogram: nugget=%.4f  psill=%.4f  range=%.1fm",
                  vgm_fit$psill[1], vgm_fit$psill[2], vgm_fit$range[2]))

  message("  Running universal kriging...")
  krige_out <- gstat::krige(
    formula   = y ~ water,
    locations = obs_df,
    newdata   = pred_df,
    model     = vgm_fit
  )

  mu_t  <- krige_out$var1.pred
  se_t  <- sqrt(pmax(krige_out$var1.var, 0))
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t
  ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }
  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  mean_SE=%.4f",
                  rmse, r2, coverage(which(ns_t >= 1L)), mean(se_t)))

  results_kriging_albedo_150m <- list(
    run_name              = "kriging_albedo_150m",
    tags                  = list(resolution = 150L, method = "universal_kriging",
                                 covariates = "water", support = "point_centroid",
                                 vgm_model = "Exp"),
    timestamp             = Sys.time(),
    resolution_m          = 150L,
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    vgm_fit               = vgm_fit,
    vgm_emp               = vgm_emp,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_kriging_albedo_150m, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------

cat("\n=== 150m albedo comparison ===\n")
cat(sprintf("%-35s  %6s  %6s  %8s  %8s\n",
            "Method", "RMSE", "R2", "cov_obs", "mean_SE"))
cat(strrep("-", 72), "\n")

for (nm in c("results_water_rsr_rho_cv_150m",
             "results_albedo_gA_centroid_150m",
             "results_kriging_albedo_150m")) {
  r <- tryCatch({
    e <- new.env()
    data(list = nm, package = "goebel2026", envir = e)
    e[[nm]]
  }, error = function(e) NULL)
  if (!is.null(r))
    cat(sprintf("%-35s  %6.4f  %6.4f  %8.3f  %8.4f\n",
                nm, r$rmse, r$r2, r$coverage_95_obs,
                mean(r$posterior_se, na.rm = TRUE)))
}

message("\nrun_150m_albedo.R complete.")
