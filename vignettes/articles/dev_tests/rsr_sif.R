## =============================================================================
## compare_rsr_sif.R
##
## Same as compare_rsr.R but using SIF_757nm as the response.
## No ground truth available -- metrics and residual plots are omitted.
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)

set.seed(42)

# =============================================================================
# 1. Data
# =============================================================================

y_alb <- goebel2026::soundings_augmented$SIF_757nm

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)

p <- ncol(A)
W <- goebel2026::make_W_matrix(goebel2026::target_grid)

Rinv <- Matrix::diag(1/soundings_augmented$SIF_Uncertainty_757nm)

# =============================================================================
# 2. Covariates
# =============================================================================

water_sounding <- goebel2026::soundings_augmented$proportion_water
water_sounding[is.na(water_sounding)] <- 0
X_fixed <- cbind(intercept = 1, water = water_sounding)
q       <- ncol(X_fixed)

water_grid_p <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
X_grid <- cbind(intercept = 1, water = water_grid_p)

# =============================================================================
# 3. Intrinsic SAR prior (rho = 1)
# =============================================================================

rho_fixed <- 1.0
IminusW   <- Matrix::Diagonal(nrow(W)) - rho_fixed * W
Q_fixed   <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
Q_fixed   <- Matrix::drop0(Q_fixed)

cat(sprintf("Q fixed at rho=1: p=%d  nnz=%d\n", nrow(Q_fixed), nnzero(Q_fixed)))

# =============================================================================
# 4. Tune phi via CV
# =============================================================================

lambda_beta <- 0.01

A_aug <- as(cbind(A, X_fixed), "dgCMatrix")

Q_aug_fun <- function(theta) {
  Q_aug     <- Matrix::bdiag(Q_fixed, lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

cat("\nTuning phi via 5-fold CV...\n")
tuned <- tune_cv(
  y          = y_alb,
  A          = A_aug,
  Q_fun      = Q_aug_fun,
  theta_init = numeric(0),
  #X_cov      = X_fixed,
  R_inv = Rinv,
  k          = 5L,
  verbose    = TRUE
)

phi_aug <- tuned$phi
cat(sprintf("\nTuned phi = %.4f\n", phi_aug))

# =============================================================================
# 5. Fit both models
# =============================================================================

Q_aug <- Matrix::bdiag(
  Q_fixed,
  lambda_beta / phi_aug * Matrix::Diagonal(q)
)

# ---- Unconstrained ----------------------------------------------------------
fit_base <- fastblm::fit_fastblm(
  y      = y_alb,
  A      = A_aug,
  Q      = Q_aug,
  phi    = phi_aug,
  R_inv = Rinv,
  solver = "cholesky"
)

beta_hat <- fit_base$posterior_mean[p + seq_len(q)]
cat(sprintf("\nUnconstrained beta: intercept=%.4f  water=%.4f\n",
            beta_hat[1], beta_hat[2]))

# ---- RSR-constrained --------------------------------------------------------
C_spatial <- as.matrix(t(X_fixed) %*% A)
C_aug     <- cbind(C_spatial, matrix(0, q, q))

fit_rsr <- fastblm::constrain(fit_base, C_aug)

beta_hat_rsr <- fit_rsr$posterior_mean[p + seq_len(q)]
cat(sprintf("RSR beta:          intercept=%.4f  water=%.4f\n",
            beta_hat_rsr[1], beta_hat_rsr[2]))
cat(sprintf("RSR ||C r||_inf = %.2e\n", max(abs(C_spatial %*% fit_rsr$posterior_mean[seq_len(p)]))))

# ---- Grid-level predictions -------------------------------------------------
A_pred   <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred     <- as.vector(A_pred %*% fit_base$posterior_mean)
pred_rsr <- as.vector(A_pred %*% fit_rsr$posterior_mean)

# =============================================================================
# 6. Plots
# =============================================================================

observed_cells <- colSums(A) > 0

plot_sf          <- target_proj
plot_sf$water    <- ifelse(observed_cells, water_grid_p, NA)
plot_sf$pred     <- ifelse(observed_cells, pred,         NA)
plot_sf$pred_rsr <- ifelse(observed_cells, pred_rsr,     NA)

sf_plot <- function(data, fill, title,
                    lo = "#2c7bb6", mi = "#ffffbf", hi = "#d7191c",
                    midpt = NULL, lims = NULL) {
  if (is.null(midpt)) midpt <- mean(data[[fill]], na.rm = TRUE)
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill]]), colour = NA) +
    scale_fill_gradient2(low = lo, mid = mi, high = hi,
                         midpoint = midpt, limits = lims,
                         na.value = "grey90", name = NULL) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          legend.key.height = unit(0.5, "cm"))
}

lim_sif <- range(c(pred, pred_rsr), na.rm = TRUE)

p1 <- sf_plot(plot_sf, "pred",     "Unconstrained",
              lims = lim_sif, midpt = mean(lim_sif))
p2 <- sf_plot(plot_sf, "pred_rsr", "RSR",
              lims = lim_sif, midpt = mean(lim_sif))
p3 <- sf_plot(plot_sf, "water",    "Water fraction",
              lo = "#f7fbff", mi = "#6baed6", hi = "#08306b", midpt = 0.25)

print(
  (p1 | p2 | p3) +
    plot_annotation(title = sprintf("Posterior mean SIF 757nm  (phi=%.4f)", phi_aug))
)

# Difference: RSR minus unconstrained
plot_sf$diff <- ifelse(observed_cells, pred_rsr - pred, NA)
lim_diff     <- max(abs(plot_sf$diff), na.rm = TRUE) * c(-1, 1)

print(
  sf_plot(plot_sf, "diff", "RSR minus unconstrained",
          lims = lim_diff, midpt = 0) +
    labs(title = "Difference: RSR − unconstrained")
)
