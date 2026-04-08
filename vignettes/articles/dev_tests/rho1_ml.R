## =============================================================================
## ml_water_rho1.R
##
## Same as ml_water.R but with rho fixed at 1 (intrinsic SAR prior) and
## phi tuned via ML. Replicates the paper setup more closely.
##
## rho=1 makes Q singular (intrinsic prior) -- K_aug is still invertible
## because A'A regularises the system in observed directions.
##
## Saves two ablations:
##   "ml_water_rho1"      -- joint fit, no RSR
##   "ml_water_rho1_rsr"  -- joint fit, RSR constraint
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)

#source("tune_ml.R")
#source("diagnostics.R")

set.seed(42)

# =============================================================================
# 1. Data setup
# =============================================================================

noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)

p <- ncol(A)

W <- goebel2026::make_W_matrix(goebel2026::target_grid)

truth <- goebel2026::target_grid$mean_albedo

# =============================================================================
# 2. Covariates
# =============================================================================

water_sounding <- goebel2026::soundings_augmented$proportion_water
water_sounding[is.na(water_sounding)] <- 0
X_fixed <- cbind(intercept = 1, water = water_sounding)   # n x 2
q       <- ncol(X_fixed)

water_grid_p <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
X_grid <- cbind(intercept = 1, water = water_grid_p)       # p x 2

# =============================================================================
# 3. Fix rho = 1, build Q (intrinsic SAR)
# =============================================================================
# rho=1 gives the intrinsic SAR prior used in the paper.
# Q = (I - W)'(I - W) is positive semidefinite (singular) but K_aug is
# still invertible because A'A is positive definite in the observed subspace.

rho_fixed <- 1.0
IminusW   <- Matrix::Diagonal(nrow(W)) - rho_fixed * W
Q_fixed   <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
Q_fixed   <- Matrix::drop0(Q_fixed)

cat(sprintf("Q fixed at rho=1: p=%d  nnz=%d\n", nrow(Q_fixed), nnzero(Q_fixed)))

# logdet of intrinsic Q is -Inf (singular), so we profile phi analytically
# using only the quadratic form -- equivalent to treating logdet_Q as a
# constant offset that doesn't affect the phi optimum.
# For the ll we use logdet_Q = 0 as a convention (improper prior).
Q_fun_rho1 <- function(theta) {
  list(Q = Q_fixed, log_det_Q = 0)   # improper prior: logdet undefined, use 0
}

# =============================================================================
# 4. Profile phi via ML with fixed rho=1
# =============================================================================
lambda_beta = 0.01

A_aug <- as(cbind(A, X_fixed), "dgCMatrix")

Q_aug_fun <- function(theta) {
  Q_aug <- Matrix::bdiag(Q_fixed, lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)   # = 2 * log(0.01) for intercept + water
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

tuned_ml_rho1 <- tune_ml(
  y          = y_alb,
  A          = A_aug,       # augmented
  Q_fun      = Q_aug_fun,   # augmented
  X_fixed    = NULL,        # no separate fixed effects -- already in A_aug
  theta_init = numeric(0),
  verbose    = TRUE
)

# =============================================================================
# 5. Build augmented system and fit
# =============================================================================

phi_aug <- tuned_ml_rho1$phi

A_aug <- as(cbind(A, X_fixed), "dgCMatrix")

# Small ridge prior on beta for RSR conditioning (see ml_water.R)
#lambda_beta <- 0.01#0.01#0.01   # desired prior precision on beta, same as old code
phi_aug     <- tuned_ml_rho1$phi

Q_aug <- Matrix::bdiag(
  Q_fixed,
  lambda_beta/phi_aug * Matrix::Diagonal(q)   # drop the phi_aug factor
)



fit_aug <- fastblm::fit_fastblm(
  y      = y_alb,
  A      = A_aug,
  Q      = Q_aug,
  phi    = phi_aug,
  solver = "cholesky"
)

r_hat    <- fit_aug$posterior_mean[seq_len(p)]
beta_hat <- fit_aug$posterior_mean[p + seq_len(q)]

cat(sprintf("\nbeta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat[1], beta_hat[2]))

se_aug  <- fastblm::posterior_se(fit_aug)
se_r    <- se_aug[seq_len(p)]
se_beta <- se_aug[p + seq_len(q)]

cat(sprintf("se_beta:  intercept=%.4f  water=%.4f\n",
            se_beta[1], se_beta[2]))

A_pred  <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
se_pred <- fastblm::posterior_se(fit_aug, A_new = A_pred)

# =============================================================================
# 6. RSR constraint
# =============================================================================

C_spatial <- as.matrix(t(X_fixed) %*% A)
C_aug     <- cbind(C_spatial, matrix(0, q, q))

fit_aug_rsr <- fastblm::constrain(fit_aug, C_aug)

r_hat_rsr    <- fit_aug_rsr$posterior_mean[seq_len(p)]
beta_hat_rsr <- fit_aug_rsr$posterior_mean[p + seq_len(q)]

cat(sprintf("\nRSR beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_rsr[1], beta_hat_rsr[2]))
cat(sprintf("RSR ||C r||_inf = %.2e\n",
            max(abs(C_spatial %*% r_hat_rsr))))

# Sanity check on C Sigma C' conditioning
SigmaCt_rsr  <- fit_aug_rsr$constraint$SigmaCt
CSigmaCt_rsr <- as.matrix(C_aug) %*% as.matrix(SigmaCt_rsr)
cat(sprintf("C Sigma C' eigenvalues: %s\n",
            paste(formatC(eigen(CSigmaCt_rsr)$values, format="e", digits=3),
                  collapse=", ")))

se_pred_rsr <- fastblm::posterior_se(fit_aug_rsr, A_new = A_pred)

# =============================================================================
# 7. Save ablations
# =============================================================================

save_ablation(
  model = list(
    fit          = fit_aug,
    se           = se_pred,
    tuned        = tuned_ml_rho1,
    X_grid       = X_grid,
    p_spatial    = p,
    colour_var   = water_grid_p,
    colour_label = "Prop. water"
  ),
  run_name = "ml_water_rho1",
  tags = list(
    tuning     = "ml",
    covariates = "intercept+water",
    constraint = "none",
    rho        = "fixed_1",
    folds      = NA
  ),
  truth   = truth,
  grid_sf = target_proj,
  A       = A
)

save_ablation(
  model = list(
    fit          = fit_aug_rsr,
    se           = se_pred_rsr,
    tuned        = tuned_ml_rho1,
    X_grid       = X_grid,
    p_spatial    = p,
    colour_var   = water_grid_p,
    colour_label = "Prop. water",
    C_check      = C_spatial
  ),
  run_name = "ml_water_rho1_rsr",
  tags = list(
    tuning     = "ml",
    covariates = "intercept+water",
    constraint = "RSR",
    rho        = "fixed_1",
    folds      = NA
  ),
  truth   = truth,
  grid_sf = target_proj,
  A       = A
)

# =============================================================================
# 8. Compare rho=1 runs
# =============================================================================

diag_list <- load_ablations(filter = list(rho = "fixed_1"))
compare_ablations(diag_list, tag = "ml_water_rho1")
