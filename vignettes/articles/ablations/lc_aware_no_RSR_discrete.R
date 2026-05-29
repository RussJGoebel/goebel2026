# =============================================================================
# Experiment: LC-aware W (binary water, alpha tuned), no RSR
# -- semi-synthetic albedo
#
# Edge weight: 1 - alpha * 1(x_i > 0.5 XOR x_j > 0.5)
# i.e. alpha penalises edges that cross the water/land boundary.
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

albedo_grid    <- goebel2026::target_grid$mean_albedo
noise_sd       <- 0.05 * sd(albedo_grid, na.rm = TRUE)
y_albedo       <- as.numeric(A %*% albedo_grid) + rnorm(nrow(A), sd = noise_sd)

cat(sprintf("Albedo noise SD: %.6f\n", noise_sd))

# -----------------------------------------------------------------------------
# LC-aware W (binary)
# -----------------------------------------------------------------------------
W_std      <- goebel2026::make_W_matrix(goebel2026::target_grid)
is_water   <- water_grid_p > 0.5   # length-p logical vector

make_W_lc_binary <- function(alpha, is_water, W_template) {
  W_trip    <- Matrix::summary(W_template)
  # 1 if edge crosses the water/land boundary, 0 if same class
  cross     <- as.numeric(is_water[W_trip$i] != is_water[W_trip$j])
  raw_w     <- pmax(1 - alpha * cross, 0)
  row_s_tab <- tapply(raw_w, W_trip$i, sum)
  row_s_vec <- numeric(nrow(W_template))
  row_s_vec[as.integer(names(row_s_tab))] <- as.numeric(row_s_tab)
  norm_w    <- raw_w / pmax(row_s_vec, .Machine$double.eps)[W_trip$i]
  Matrix::sparseMatrix(i = W_trip$i, j = W_trip$j,
                       x = as.numeric(norm_w), dims = dim(W_template))
}

make_Q_spatial <- function(alpha) {
  IminusW <- Matrix::Diagonal(p) - make_W_lc_binary(alpha, is_water, W_std)
  Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))
}

# -----------------------------------------------------------------------------
# Augmented design matrix and Q_fun
# -----------------------------------------------------------------------------
lambda_beta <- 0.01
A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

Q_fun <- function(theta) {
  Q_aug     <- Matrix::bdiag(make_Q_spatial(theta[["alpha"]]),
                             lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

# -----------------------------------------------------------------------------
# CV
# -----------------------------------------------------------------------------
cat("Tuning LC-aware W binary (alpha in [0,1]), no RSR ...\n")
t_cv_start <- proc.time()
tuned <- tune_cv(
  y          = y_albedo,
  A          = A_aug,
  Q_fun      = Q_fun,
  theta_init = c(alpha = 0.5),
  lower      = c(alpha = 0.0),
  upper      = c(alpha = 1.0),
  k          = 10L,
  constraint = NULL,
  verbose    = TRUE
)
t_cv <- proc.time() - t_cv_start
cat(sprintf("alpha = %.4f   phi = %.4f   CV-MSE = %.6f\n",
            tuned$theta[["alpha"]], tuned$phi, tuned$value))
cat(sprintf("CV time: %.1fs\n", t_cv["elapsed"]))

# -----------------------------------------------------------------------------
# Final fit at tuned parameters
# -----------------------------------------------------------------------------
t_fit_start <- proc.time()
alpha_hat <- tuned$theta[["alpha"]]
phi_hat   <- tuned$phi

Q_sp  <- make_Q_spatial(alpha_hat)
Q_aug <- Matrix::bdiag(Q_sp, lambda_beta / phi_hat * Matrix::Diagonal(q))

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
# Save
# -----------------------------------------------------------------------------
t_total <- proc.time() - t_start
timings <- list(cv = t_cv, fit = t_fit, total = t_total)
cat(sprintf("Total time: %.1fs\n", t_total["elapsed"]))

results_lc_binary_norsr_albedo <- list(
  tuned      = tuned,
  fit        = fit,
  alpha_hat  = alpha_hat,
  phi_hat    = phi_hat,
  beta_hat   = beta_hat,
  noise_sd   = noise_sd,
  y_albedo   = y_albedo,
  timings    = timings,
  p          = p,
  q          = q
)

usethis::use_data(results_lc_binary_norsr_albedo, overwrite = TRUE)
cat("Saved via use_data\n")
