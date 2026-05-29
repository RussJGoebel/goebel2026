# =============================================================================
# Experiment: land-fraction weighted Q (no RSR) -- semi-synthetic albedo
#
# Prior on spatial field r is attenuated by land fraction:
#   Q_pinned = D_land %*% Q_sp %*% D_land
# where D_land = diag(1 - water_grid_p).
# Pure water cells contribute nothing to r; mixtures interpolate.
# =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)

set.seed(42)
t_start <- proc.time()

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)
p <- ncol(A)

water_sounding <- goebel2026::soundings_augmented$proportion_water
X_fixed        <- cbind(intercept = 1, water = water_sounding)
q              <- ncol(X_fixed)

water_grid_p   <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
land_frac      <- 1 - water_grid_p

albedo_grid    <- goebel2026::target_grid$mean_albedo
noise_sd       <- 0.05 * sd(albedo_grid, na.rm = TRUE)
y_albedo       <- as.numeric(A %*% albedo_grid) + rnorm(nrow(A), sd = noise_sd)

cat(sprintf("Albedo noise SD: %.6f\n", noise_sd))

# -----------------------------------------------------------------------------
# W and Q construction
# -----------------------------------------------------------------------------
W_std <- goebel2026::make_W_matrix(goebel2026::target_grid)

# Precision bump per unit water fraction. Controls how strongly water
# cells are pulled toward zero. Tunable but not CV'd here.
water_kappa <- 100

make_Q_spatial <- function(alpha = 0) {
  IminusW <- Matrix::Diagonal(p) - W_std
  Q_sp    <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))

  # Add precision proportional to water fraction on the diagonal only.
  # Pure water cells get +water_kappa precision (strongly pinned to 0);
  # pure land cells are unaffected; mixtures interpolate.
  # Spatial connectivity (off-diagonal) is left intact.
  Matrix::drop0(Q_sp + Matrix::Diagonal(x = water_kappa * water_grid_p))
}

# -----------------------------------------------------------------------------
# Augmented design matrix and Q_fun
# -----------------------------------------------------------------------------
lambda_beta <- 0.01
A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

# No theta to tune (standard W + fixed land-fraction weighting);
# phi only via CV.
Q_fun <- function(theta) {
  Q_aug     <- Matrix::bdiag(make_Q_spatial(),
                             lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

# -----------------------------------------------------------------------------
# CV
# -----------------------------------------------------------------------------
cat("Tuning phi (land-fraction weighted Q, no RSR) ...\n")
t_cv_start <- proc.time()
tuned <- tune_cv(
  y          = y_albedo,
  A          = A_aug,
  Q_fun      = Q_fun,
  theta_init = numeric(0),
  k          = 10L,
  constraint = NULL,
  verbose    = TRUE
)
t_cv <- proc.time() - t_cv_start
cat(sprintf("phi = %.4f   CV-MSE = %.6f\n", tuned$phi, tuned$value))
cat(sprintf("CV time: %.1fs\n", t_cv["elapsed"]))

# -----------------------------------------------------------------------------
# Final fit
# -----------------------------------------------------------------------------
t_fit_start <- proc.time()
phi_hat <- tuned$phi

Q_sp_pinned <- make_Q_spatial()
Q_aug       <- Matrix::bdiag(Q_sp_pinned,
                             lambda_beta / phi_hat * Matrix::Diagonal(q))

fit <- fastblm::fit_fastblm(
  y      = y_albedo,
  A      = A_aug,
  Q      = Q_aug,
  phi    = phi_hat,
  solver = "cholesky"
)
t_fit <- proc.time() - t_fit_start

beta_hat <- fit$posterior_mean[p + seq_len(q)]
cat(sprintf("beta: intercept = %.4f   water = %.4f\n",
            beta_hat[1], beta_hat[2]))
cat(sprintf("Fit time: %.1fs\n", t_fit["elapsed"]))

# -----------------------------------------------------------------------------
# Posterior mean and SE on grid
# -----------------------------------------------------------------------------
t_pred_start <- proc.time()
X_grid    <- cbind(intercept = 1, water = water_grid_p)
A_pred    <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred_grid <- as.numeric(A_pred %*% fit$posterior_mean)
se_grid   <- fastblm::posterior_se(fit, A_new = A_pred)
t_pred    <- proc.time() - t_pred_start
cat(sprintf("Prediction time: %.1fs\n", t_pred["elapsed"]))

# -----------------------------------------------------------------------------
# Summary statistics on observed cells
# -----------------------------------------------------------------------------
observed   <- which(Matrix::colSums(A) > 0)
true_grid  <- albedo_grid

resid      <- pred_grid[observed] - true_grid[observed]
ss_res     <- sum(resid^2)
ss_tot     <- sum((true_grid[observed] - mean(true_grid[observed]))^2)
r2         <- 1 - ss_res / ss_tot
rmse       <- sqrt(mean(resid^2))

# 95% credible interval coverage
lower_95   <- pred_grid[observed] - 1.96 * se_grid[observed]
upper_95   <- pred_grid[observed] + 1.96 * se_grid[observed]
coverage   <- mean(true_grid[observed] >= lower_95 &
                     true_grid[observed] <= upper_95)

cat(sprintf("\n--- Summary (observed cells, n=%d) ---\n", length(observed)))
cat(sprintf("RMSE     : %.6f\n", rmse))
cat(sprintf("R2       : %.4f\n", r2))
cat(sprintf("Coverage : %.4f (nominal 0.95)\n", coverage))

stats <- list(rmse = rmse, r2 = r2, coverage = coverage, n_observed = length(observed))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
t_total <- proc.time() - t_start
timings <- list(cv = t_cv, fit = t_fit, pred = t_pred, total = t_total)
cat(sprintf("Total time: %.1fs\n", t_total["elapsed"]))

results_landfrac_norsr_albedo <- list(
  tuned      = tuned,
  fit        = fit,
  phi_hat    = phi_hat,
  beta_hat   = beta_hat,
  pred_grid  = pred_grid,
  se_grid    = se_grid,
  stats      = stats,
  noise_sd   = noise_sd,
  y_albedo   = y_albedo,
  timings    = timings,
  p          = p,
  q          = q
)

usethis::use_data(results_landfrac_norsr_albedo, overwrite = TRUE)
cat("Saved via use_data\n")
