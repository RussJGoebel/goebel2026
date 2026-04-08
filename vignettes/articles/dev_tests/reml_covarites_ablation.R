## =============================================================================
## ml_water.R
##
## Marginal ML tuning with intercept + water as covariates, then joint fit
## of (spatial field, beta) via the augmented system [A | X_fixed].
##
## Model:  y = A r + X_fixed beta + eps
##         r    ~ N(0, sigma2b * Q^{-1})   SAR prior
##         beta ~ flat (improper, zero precision block)
##
## Hyperparameters (rho, phi) estimated by tune_ml() -- exact sparse Cholesky
## on the augmented K_aug, no Lanczos, no projection.
##
## Posterior SEs for diagnostics are computed for the predicted field
##   yhat_i = r_i + beta_0 + beta_water * water_grid_i
## via posterior_se(fit, A_new = [I_p | X_grid]), which accounts for joint
## uncertainty in both r and beta.
##
## Saves two ablations:
##   "ml_water"      -- joint fit, no RSR constraint
##   "ml_water_rsr"  -- joint fit, RSR constraint on spatial block
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)

#source("tune_ml.R")       # tune_ml() -- add to fastblm package when ready
source("diagnostics.R")   # blm_diagnostics, save_ablation, etc.

set.seed(42)

# =============================================================================
# 1. Data setup
# =============================================================================

noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)

p <- ncol(A)   # number of spatial grid cells

W <- goebel2026::make_W_matrix(goebel2026::target_grid)

Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  Q       <- Matrix::drop0(Q)
  # Exact logdet via sparse Cholesky -- avoids Lanczos bias entirely
  CQ      <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
  ld_Q    <- as.numeric(
    Matrix::determinant(CQ, logarithm = TRUE, sqrt = TRUE)$modulus) * 2
  list(Q = Q, log_det_Q = ld_Q)
}

truth <- goebel2026::target_grid$mean_albedo

# =============================================================================
# 2. Covariates
# =============================================================================
# X_fixed: covariates at the sounding level (n x 2) -- used for tuning/fitting.
# X_grid:  covariates at the grid level (p x 2) -- used for prediction SEs.

water_sounding <- goebel2026::soundings_augmented$proportion_water
water_sounding[is.na(water_sounding)] <- 0
X_fixed <- cbind(intercept = 1, water = water_sounding)   # n x 2
q       <- ncol(X_fixed)

water_grid_p <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
X_grid <- cbind(intercept = 1, water = water_grid_p)       # p x 2

# =============================================================================
# 3. Tune rho and phi via marginal ML
# =============================================================================
# tune_ml augments A with X_fixed internally, builds Q_aug = bdiag(Q, 0),
# and uses sparse Cholesky of K_aug at every evaluation.
# No Lanczos -- logdet exact from Cholesky. No projection -- ML not REML.

tuned_ml_water <- tune_ml(
  y          = y_alb,
  A          = A,
  Q_fun      = Q_fun,
  X_fixed    = X_fixed,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 0.99,
  verbose    = TRUE
)

cat(sprintf("\nML+water optimum: rho=%.4f  phi=%.4f  sigma2e=%.6g\n",
            tuned_ml_water$theta[["rho"]],
            tuned_ml_water$phi,
            tuned_ml_water$sigma2e))

# =============================================================================
# 4. Build augmented system for the final fit
# =============================================================================

rho_opt <- tuned_ml_water$theta[["rho"]]
Q_opt   <- Matrix::forceSymmetric(
  Matrix::crossprod(Matrix::Diagonal(nrow(W)) - rho_opt * W)
)
phi_aug <- tuned_ml_water$phi

A_aug <- as(cbind(A, X_fixed), "dgCMatrix")   # n x (p+q)

# Small ridge prior on beta: necessary for RSR to be well-conditioned.
# A flat prior (zero precision) inflates the marginal posterior covariance
# of r in the beta-correlated directions by ~5000x, making C Sigma C'
# nearly singular and the RSR correction numerically catastrophic.
# lambda_beta = 1e-6 is weak enough to leave beta estimates unchanged
# but bounds Sigma_rr so that constrain() works correctly -- matching
# what the original fitted_spatial_model code did via compute_precision().
lambda_beta <- 1e-6

Q_aug <- Matrix::forceSymmetric(Matrix::bdiag(
  Q_opt,
  lambda_beta * Matrix::Diagonal(q)
))

# =============================================================================
# 5. Joint fit of (r, beta) -- exact Cholesky solve
# =============================================================================
# solver="cholesky" caches the Cholesky factor, enabling exact posterior SEs
# via the Cholesky path in posterior_se() -- no Hutchinson approximation.

fit_aug <- fastblm::fit_fastblm(
  y      = y_alb,
  A      = A_aug,
  Q      = Q_aug,
  phi    = phi_aug,
  solver = "cholesky"
)

# Split posterior mean into spatial field and beta
r_hat    <- fit_aug$posterior_mean[seq_len(p)]
beta_hat <- fit_aug$posterior_mean[p + seq_len(q)]

cat(sprintf("\nbeta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat[1], beta_hat[2]))

# --- Component SEs (for reporting) ------------------------------------------
se_aug  <- fastblm::posterior_se(fit_aug)        # SE of [r; beta]
se_r    <- se_aug[seq_len(p)]
se_beta <- se_aug[p + seq_len(q)]

cat(sprintf("se_beta:  intercept=%.4f  water=%.4f\n",
            se_beta[1], se_beta[2]))

# --- Prediction SEs (for diagnostics / coverage) ----------------------------
# SE of yhat_i = r_i + x_grid_i' beta, accounting for joint uncertainty in
# both r and beta. A_pred = [I_p | X_grid] so each row gives one prediction.
A_pred  <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")   # p x (p+q)
se_pred <- fastblm::posterior_se(fit_aug, A_new = A_pred)         # p-vector

# =============================================================================
# 6. RSR constraint on the spatial block
# =============================================================================
# C_spatial = t(X_fixed) %*% A  (2 x p): forces r orthogonal to intercept+water
# as seen through A. Pad with zeros for the beta columns.

C_spatial <- as.matrix(t(X_fixed) %*% A)          # 2 x p
C_aug     <- cbind(C_spatial, matrix(0, q, q))     # 2 x (p+q)

fit_aug_rsr <- fastblm::constrain(fit_aug, C_aug)

r_hat_rsr    <- fit_aug_rsr$posterior_mean[seq_len(p)]
beta_hat_rsr <- fit_aug_rsr$posterior_mean[p + seq_len(q)]

cat(sprintf("\nRSR beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_rsr[1], beta_hat_rsr[2]))
cat(sprintf("RSR constraint residual ||C r||_inf = %.2e\n",
            max(abs(C_spatial %*% r_hat_rsr))))

# RSR prediction SEs -- posterior_se applies Schur correction automatically
se_pred_rsr <- fastblm::posterior_se(fit_aug_rsr, A_new = A_pred)   # p-vector

# =============================================================================
# 7. Save both ablations
# =============================================================================

save_ablation(
  model = list(
    fit          = fit_aug,
    se           = se_pred,         # SE of predicted field r + X_grid beta
    tuned        = tuned_ml_water,
    X_grid       = X_grid,          # needed by blm_diagnostics for mu_pred
    p_spatial    = p,
    colour_var   = water_grid_p,
    colour_label = "Prop. water"
  ),
  run_name = "ml_water",
  tags = list(
    tuning     = "ml",
    covariates = "intercept+water",
    constraint = "none",
    folds      = NA
  ),
  truth   = truth,
  grid_sf = target_proj,
  A       = A
)

save_ablation(
  model = list(
    fit          = fit_aug_rsr,
    se           = se_pred_rsr,     # SE of predicted field under RSR
    tuned        = tuned_ml_water,
    X_grid       = X_grid,
    p_spatial    = p,
    colour_var   = water_grid_p,
    colour_label = "Prop. water",
    C_check      = C_spatial        # constraint residual check on r only
  ),
  run_name = "ml_water_rsr",
  tags = list(
    tuning     = "ml",
    covariates = "intercept+water",
    constraint = "RSR",
    folds      = NA
  ),
  truth   = truth,
  grid_sf = target_proj,
  A       = A
)

# =============================================================================
# 8. Compare
# =============================================================================

diag_list <- load_ablations(filter = list(tuning = "ml"))
compare_ablations(diag_list, tag = "ml_water")
