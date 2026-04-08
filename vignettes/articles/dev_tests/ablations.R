## =============================================================================
## Ablation experiments 1-3
## =============================================================================
## 1. Spatial k-means CV folds (vs. random folds)
## 2. Covariates via REML projection (X_fixed in tune_reml)
## 3. Constrained fit with covariates (RSR: C = t(X_cov) %*% A)
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(future)

set.seed(42)

# -----------------------------------------------------------------------------
# Shared setup  (same as baseline example)
# -----------------------------------------------------------------------------
noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

soundings_proj  <- spatintegrate::ensure_projected(goebel2026::soundings)
target_proj     <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- spatintegrate::compute_overlap_fractions(soundings_proj, target_proj)
A <- as(A, "dgCMatrix")

W       <- goebel2026::make_W_matrix(goebel2026::target_grid)
Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  list(Q = Matrix::drop0(Q))
}

truth   <- goebel2026::target_grid$mean_albedo
na_idx  <- is.na(truth)

# Helper: evaluate a fit + se against truth
eval_fit <- function(fit, se, label) {
  rmse <- sqrt(mean((fit$posterior_mean[!na_idx] - truth[!na_idx])^2))
  lower <- fit$posterior_mean - 1.96 * se
  upper <- fit$posterior_mean + 1.96 * se
  cov  <- mean(truth[!na_idx] >= lower[!na_idx] &
                 truth[!na_idx] <= upper[!na_idx])
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  rho:      %.4f\n", fit$phi))   # phi from fit; rho from tuned
  cat(sprintf("  phi:      %.4f\n", fit$phi))
  cat(sprintf("  sigma2e:  %.4f\n", fit$sigma2e))
  cat(sprintf("  RMSE:     %.4f\n", rmse))
  cat(sprintf("  Coverage: %.4f\n", cov))
  invisible(list(rmse = rmse, coverage = cov))
}

# =============================================================================
# Ablation 1: Spatial k-means CV folds
# =============================================================================
# Assign soundings to folds by clustering their centroids, so that
# held-out folds are spatially separated from training — a harder and more
# realistic CV that penalises spatial overfitting.

make_spatial_folds <- function(sf_obj, k, seed = 42) {
  cents  <- sf::st_centroid(sf_obj)
  coords <- sf::st_coordinates(cents)   # n x 2 matrix (X, Y)
  set.seed(seed)
  km     <- kmeans(coords, centers = k, nstart = 25, iter.max = 100)
  as.integer(km$cluster)
}

k <- 10L
spatial_folds <- make_spatial_folds(soundings_proj, k = k)

cat("Fold sizes (spatial k-means):\n")
print(table(spatial_folds))

future::plan(future::multisession())

tuned_cv_spatial <- fastblm::tune_cv(
  y          = y_alb,
  A          = A,
  Q_fun      = Q_fun,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 0.999,
  folds      = spatial_folds,   # <-- pre-computed spatial folds
  solver     = "cholesky",
  parallel   = TRUE,
  verbose    = TRUE
)

future::plan(future::sequential())

# Refit at tuned hyperparameters using the tuned rho
rho_cv_sp <- tuned_cv_spatial$theta[["rho"]]
IminusW_cv_sp <- Matrix::Diagonal(nrow(W)) - rho_cv_sp * W
Q_cv_sp <- Matrix::forceSymmetric(Matrix::crossprod(IminusW_cv_sp))

fit_cv_spatial <- fastblm::fit_fastblm(
  y_alb, A, Q_cv_sp, tuned_cv_spatial$phi, solver = "cholesky"
)
se_cv_spatial <- fastblm::posterior_se(fit_cv_spatial)

cat(sprintf("\nSpatial CV: rho=%.4f  phi=%.4f\n",
            tuned_cv_spatial$theta[["rho"]], tuned_cv_spatial$phi))
eval_fit(fit_cv_spatial, se_cv_spatial, "Ablation 1: Spatial k-means CV")

# =============================================================================
# Ablation 2: Covariates via REML projection  (X_fixed in tune_reml)
# =============================================================================
# Use water fraction as a covariate. REML projects it out of y and A before
# optimising, giving variance component estimates free of covariate confounding.
# The posterior mean is then recovered by adding back the fixed-effect fit.

# Build covariate matrix at the sounding level: X_fixed is n x q
# water_frac upscaled to soundings via A (same overlap matrix)
water_grid   <- goebel2026::target_grid$prop_water          # p-vector on grid
water_grid[is.na(water_grid)] <- 0
water_sounding <- as.numeric(A %*% water_grid)              # n-vector

X_fixed <- cbind(intercept = 1, water = water_sounding)    # n x 2

tuned_reml_cov <- fastblm::tune_reml(
  y             = y_alb,
  A             = A,
  Q_fun         = Q_fun,
  X_fixed       = X_fixed,
  theta_init    = c(rho = 0.9),
  lower         = 0.01,
  upper         = 0.999,
  solver        = "cholesky",
  logdet_method = "cholesky",
  verbose       = TRUE
)

rho_reml_cov <- tuned_reml_cov$theta[["rho"]]
IminusW_reml_cov <- Matrix::Diagonal(nrow(W)) - rho_reml_cov * W
Q_reml_cov <- Matrix::forceSymmetric(Matrix::crossprod(IminusW_reml_cov))

# For the final fit we include fixed effects explicitly by appending columns.
# Strategy: augment A with X_fixed columns and a block-diagonal Q that
# puts a flat (large-variance) prior on the fixed-effect coefficients,
# OR simply fit with the projected y/A directly and recover beta separately.
#
# Here we use the simpler approach: project out fixed effects from y and A,
# fit the spatial field on residuals, then report residual-space metrics.
proj <- fastblm:::reml_project(y_alb, A, X_fixed)   # internal helper

fit_reml_cov <- fastblm::fit_fastblm(
  proj$y, proj$A, Q_reml_cov, tuned_reml_cov$phi, solver = "cholesky"
)
se_reml_cov <- fastblm::posterior_se(fit_reml_cov)

cat(sprintf("\nREML+covariates: rho=%.4f  phi=%.4f\n",
            tuned_reml_cov$theta[["rho"]], tuned_reml_cov$phi))
# Note: RMSE/coverage are in the projected (residual) space here, so values
# are not directly comparable to the baseline; see note below.
rmse_cov <- sqrt(mean((fit_reml_cov$posterior_mean[!na_idx] -
                         truth[!na_idx])^2))
cat(sprintf("  RMSE (residual space): %.4f\n", rmse_cov))
cat(sprintf("  sigma2e: %.4f  sigma2b: %.4f\n",
            fit_reml_cov$sigma2e, fit_reml_cov$sigma2b))

# =============================================================================
# Ablation 3: Constrained fit with covariates  (RSR)
# =============================================================================
# Restricted Spatial Regression: enforce X^(o)' A gamma = 0, i.e. the spatial
# field is orthogonal to the upscaled covariates. This is the paper's
# Section 2.5.1 approach.
#
# C = t(X_fixed) %*% A   (q x p)
# constrain(fit, C) applies the Schur-complement correction in-place.

# Step 1: tune hyperparameters with the constraint enforced on each CV fold
# using X_cov so the fold-level C is recomputed per fold (correct approach).

future::plan(future::multisession())

tuned_cv_rsr <- fastblm::tune_cv(
  y          = y_alb,
  A          = A,
  Q_fun      = Q_fun,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 0.999,
  folds      = spatial_folds,   # reuse spatial folds for fair comparison
  X_cov      = X_fixed,         # per-fold C = t(X_fixed[train,]) %*% A_train
  solver     = "cholesky",
  parallel   = TRUE,
  verbose    = TRUE
)

future::plan(future::sequential())

# Step 2: fit on full data and apply global constraint
rho_rsr <- tuned_cv_rsr$theta[["rho"]]
IminusW_rsr <- Matrix::Diagonal(nrow(W)) - rho_rsr * W
Q_rsr <- Matrix::forceSymmetric(Matrix::crossprod(IminusW_rsr))

fit_rsr_unconstrained <- fastblm::fit_fastblm(
  y_alb, A, Q_rsr, tuned_cv_rsr$phi, solver = "cholesky"
)

# Global constraint matrix C = t(X_fixed) %*% A  (q x p)
C_global <- t(X_fixed) %*% A

fit_rsr <- fastblm::constrain(fit_rsr_unconstrained, C_global)
se_rsr  <- fastblm::posterior_se(fit_rsr)

cat(sprintf("\nRSR (constrained + covariates): rho=%.4f  phi=%.4f\n",
            tuned_cv_rsr$theta[["rho"]], tuned_cv_rsr$phi))
eval_fit(fit_rsr, se_rsr, "Ablation 3: RSR constrained")

# =============================================================================
# Summary table
# =============================================================================
results <- data.frame(
  model    = c("Baseline CV (random folds)",
               "Ablation 1: Spatial k-means CV",
               "Ablation 2: REML + covariates (residual space)",
               "Ablation 3: RSR constrained"),
  rho      = c(NA,                               # fill from baseline run
               tuned_cv_spatial$theta[["rho"]],
               tuned_reml_cov$theta[["rho"]],
               tuned_cv_rsr$theta[["rho"]]),
  phi      = c(NA,
               tuned_cv_spatial$phi,
               tuned_reml_cov$phi,
               tuned_cv_rsr$phi),
  sigma2e  = c(NA,
               tuned_cv_spatial$sigma2e,
               tuned_reml_cov$sigma2e,
               tuned_cv_rsr$sigma2e)
)
print(results, digits = 4)

usethis::use_data(results,overwrite = TRUE)
