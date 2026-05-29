# =============================================================================
# Experiment 1: Landcover-aware SAR prior
#
# Modifies the W matrix so that adjacency weights between pixels are
# downweighted proportionally to the difference in water fraction:
#
#   raw_w(i,j) = 1 - alpha * |x_i - x_j|    (clamped >= 0)
#
# then row-normalised so each row sums to 1. alpha in [0,1] is tuned by CV
# jointly with phi. alpha=0 recovers the standard uniform W.
# =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)

set.seed(42)

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

truth    <- goebel2026::target_grid$mean_albedo
noise_sd <- sd(truth, na.rm = TRUE) / 20
y_alb    <- A %*% truth + rnorm(nrow(A), 0, noise_sd)

water_sounding             <- goebel2026::soundings_augmented$proportion_water
X_fixed                    <- cbind(intercept = 1, water = water_sounding)
q                          <- ncol(X_fixed)

water_grid_p               <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
X_grid                     <- cbind(intercept = 1, water = water_grid_p)

# -----------------------------------------------------------------------------
# Landcover-aware W
#
# Uses the sparsity pattern of the standard W (queen adjacency) but replaces
# the uniform 1/n_i weights with landcover-similarity weights, then
# row-normalises.
# -----------------------------------------------------------------------------
W_std <- goebel2026::make_W_matrix(goebel2026::target_grid)

make_W_lc <- function(alpha, water_prop, W_template) {
  W_trip <- Matrix::summary(W_template)   # data.frame with columns i, j, x

  xi     <- water_prop[W_trip$i]
  xj     <- water_prop[W_trip$j]
  raw_w  <- pmax(1 - alpha * abs(xi - xj), 0)

  # tapply returns a named array; as.numeric() strips that before indexing
  row_s  <- as.numeric(tapply(raw_w, W_trip$i, sum))
  norm_w <- raw_w / row_s[W_trip$i]

  Matrix::sparseMatrix(
    i    = W_trip$i,
    j    = W_trip$j,
    x    = as.numeric(norm_w),   # ensure plain numeric, not named array
    dims = dim(W_template)
  )
}

# Helper: build Q_aug given alpha
make_Q_aug <- function(alpha, lambda_beta, q, p) {
  W_lc      <- make_W_lc(alpha, water_grid_p, W_std)
  IminusW   <- Matrix::Diagonal(p) - W_lc
  Q_sp      <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  Q_sp      <- Matrix::drop0(Q_sp)
  Q_aug     <- Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))
  Q_aug
}

# -----------------------------------------------------------------------------
# Augmented design matrix
# -----------------------------------------------------------------------------
lambda_beta <- 0.01
A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

# -----------------------------------------------------------------------------
# Q_fun for tune_cv: theta = c(alpha)
# log_det_Q convention: intrinsic SAR is improper so we use 0 for the spatial
# block; only the beta ridge contributes.
# -----------------------------------------------------------------------------
Q_aug_fun <- function(theta) {
  alpha     <- if (length(theta) > 0) theta[["alpha"]] else 0.0
  Q_aug     <- make_Q_aug(alpha, lambda_beta, q, p)
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

# -----------------------------------------------------------------------------
# Tune alpha and phi jointly via 5-fold CV
# -----------------------------------------------------------------------------
tuned_lc <- tune_cv(
  y          = y_alb,
  A          = A_aug,
  Q_fun      = Q_aug_fun,
  theta_init = c(alpha = 0.5),
  lower      = c(alpha = 0.0),
  upper      = c(alpha = 1.0),
  k          = 5L,
  verbose    = TRUE
)

alpha_hat <- tuned_lc$theta[["alpha"]]
phi_aug   <- tuned_lc$phi
cat(sprintf("\nTuned alpha = %.4f   phi = %.4f\n", alpha_hat, phi_aug))

# -----------------------------------------------------------------------------
# Fit: landcover-aware model
# -----------------------------------------------------------------------------
Q_lc  <- make_Q_aug(alpha_hat, lambda_beta, q, p)

# Re-scale beta block: the fit uses Q/phi internally, but the beta block needs
# lambda_beta/phi scaling to match the original parameterisation
Q_lc_fit <- Matrix::bdiag(
  make_Q_aug(alpha_hat, lambda_beta, q, p)[seq_len(p), seq_len(p)],
  lambda_beta / phi_aug * Matrix::Diagonal(q)
)

fit_lc <- fastblm::fit_fastblm(
  y      = y_alb,
  A      = A_aug,
  Q      = Q_lc_fit,
  phi    = phi_aug,
  solver = "cholesky"
)

r_hat    <- fit_lc$posterior_mean[seq_len(p)]
beta_hat <- fit_lc$posterior_mean[p + seq_len(q)]
cat(sprintf("beta_hat: intercept=%.4f  water=%.4f\n", beta_hat[1], beta_hat[2]))

se_lc      <- fastblm::posterior_se(fit_lc)
se_r_lc    <- se_lc[seq_len(p)]
se_beta_lc <- se_lc[p + seq_len(q)]
cat(sprintf("se_beta:  intercept=%.4f  water=%.4f\n", se_beta_lc[1], se_beta_lc[2]))

A_pred     <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred_lc    <- as.vector(A_pred %*% fit_lc$posterior_mean)
se_pred_lc <- fastblm::posterior_se(fit_lc, A_new = A_pred)

# -----------------------------------------------------------------------------
# RSR constraint
# -----------------------------------------------------------------------------
C_spatial   <- as.matrix(t(X_fixed) %*% A)
C_aug       <- cbind(C_spatial, matrix(0, q, q))

fit_lc_rsr     <- fastblm::constrain(fit_lc, C_aug)
pred_lc_rsr    <- as.vector(A_pred %*% fit_lc_rsr$posterior_mean)
se_pred_lc_rsr <- fastblm::posterior_se(fit_lc_rsr, A_new = A_pred)

beta_hat_rsr <- fit_lc_rsr$posterior_mean[p + seq_len(q)]
cat(sprintf("RSR beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_rsr[1], beta_hat_rsr[2]))
cat(sprintf("RSR ||C r||_inf = %.2e\n",
            max(abs(C_spatial %*% fit_lc_rsr$posterior_mean[seq_len(p)]))))

# -----------------------------------------------------------------------------
# Baseline: standard W (alpha = 0) at same phi for fair comparison
# -----------------------------------------------------------------------------
Q_base_fit <- Matrix::bdiag(
  make_Q_aug(0.0, lambda_beta, q, p)[seq_len(p), seq_len(p)],
  lambda_beta / phi_aug * Matrix::Diagonal(q)
)

fit_base   <- fastblm::fit_fastblm(
  y = y_alb, A = A_aug, Q = Q_base_fit, phi = phi_aug, solver = "cholesky"
)
pred_base    <- as.vector(A_pred %*% fit_base$posterior_mean)
se_pred_base <- fastblm::posterior_se(fit_base, A_new = A_pred)

fit_base_rsr     <- fastblm::constrain(fit_base, C_aug)
pred_base_rsr    <- as.vector(A_pred %*% fit_base_rsr$posterior_mean)
se_pred_base_rsr <- fastblm::posterior_se(fit_base_rsr, A_new = A_pred)

# -----------------------------------------------------------------------------
# Evaluation
# -----------------------------------------------------------------------------
observed_cells <- colSums(A) > 0

eval_metrics <- function(label, pred, se, truth) {
  keep     <- observed_cells & !is.na(truth)
  resid    <- pred[keep] - truth[keep]
  rmse     <- sqrt(mean(resid^2))
  r2       <- 1 - sum(resid^2) / sum((truth[keep] - mean(truth[keep]))^2)
  coverage <- mean(truth[keep] >= pred[keep] - 1.96 * se[keep] &
                     truth[keep] <= pred[keep] + 1.96 * se[keep])
  cat(sprintf("%-40s  RMSE=%.4f  R2=%.4f  95%%cov=%.1f%%  (n=%d)\n",
              label, rmse, r2, 100 * coverage, sum(keep)))
  data.frame(Model = label, RMSE = rmse, R2 = r2, Coverage = coverage)
}

m1 <- eval_metrics("Standard W — unconstrained",
                   pred_base,    se_pred_base,    truth)
m2 <- eval_metrics("Standard W — RSR",
                   pred_base_rsr, se_pred_base_rsr, truth)
m3 <- eval_metrics(sprintf("LC-aware W (alpha=%.3f) — unconstrained", alpha_hat),
                   pred_lc,    se_pred_lc,    truth)
m4 <- eval_metrics(sprintf("LC-aware W (alpha=%.3f) — RSR", alpha_hat),
                   pred_lc_rsr, se_pred_lc_rsr, truth)

knitr::kable(rbind(m1, m2, m3, m4), digits = 4, row.names = FALSE,
             col.names = c("Model", "RMSE", "R²", "95% Coverage"))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
plot_sf <- target_proj
plot_sf$truth        <- truth
plot_sf$water        <- water_grid_p
plot_sf$pred_base    <- pred_base
plot_sf$pred_lc      <- pred_lc
plot_sf$pred_lc_rsr  <- pred_lc_rsr
plot_sf$se_base      <- se_pred_base
plot_sf$se_lc        <- se_pred_lc
plot_sf$resid_base   <- pred_base   - truth
plot_sf$resid_lc     <- pred_lc     - truth
plot_sf$resid_lc_rsr <- pred_lc_rsr - truth

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

lim_alb <- range(c(truth, pred_base, pred_lc, pred_lc_rsr), na.rm = TRUE)
lim_res <- max(abs(c(plot_sf$resid_base, plot_sf$resid_lc,
                     plot_sf$resid_lc_rsr)), na.rm = TRUE) * c(-1, 1)
lim_se  <- range(c(se_pred_base, se_pred_lc), na.rm = TRUE)

# Fitted fields
(sf_plot(plot_sf, "truth",     "True albedo",         lims = lim_alb, midpt = mean(lim_alb)) |
    sf_plot(plot_sf, "pred_base", "Standard W",          lims = lim_alb, midpt = mean(lim_alb)) |
    sf_plot(plot_sf, "pred_lc",   sprintf("LC-aware W (α=%.2f)", alpha_hat),
            lims = lim_alb, midpt = mean(lim_alb))) +
  plot_annotation(title = sprintf("Fitted fields — phi=%.4f", phi_aug))

# Residuals
(sf_plot(plot_sf, "resid_base",   "Residual — standard W",  lims = lim_res, midpt = 0) |
    sf_plot(plot_sf, "resid_lc",     "Residual — LC-aware (unconstrained)", lims = lim_res, midpt = 0) |
    sf_plot(plot_sf, "resid_lc_rsr", "Residual — LC-aware RSR", lims = lim_res, midpt = 0))

# Posterior SE comparison
(sf_plot(plot_sf, "se_base", "SE — standard W",
         lo = "#fff7ec", mi = "#fc8d59", hi = "#7f0000",
         midpt = mean(lim_se), lims = lim_se) |
    sf_plot(plot_sf, "se_lc",   sprintf("SE — LC-aware W (α=%.2f)", alpha_hat),
            lo = "#fff7ec", mi = "#fc8d59", hi = "#7f0000",
            midpt = mean(lim_se), lims = lim_se))

# Scatter: all four models
keep  <- !is.na(truth)
sc_df <- data.frame(
  truth  = rep(truth[keep], 4),
  fitted = c(pred_base[keep], pred_base_rsr[keep],
             pred_lc[keep],   pred_lc_rsr[keep]),
  model  = rep(c("Standard W", "Standard W + RSR",
                 sprintf("LC-aware (α=%.2f)", alpha_hat),
                 sprintf("LC-aware RSR (α=%.2f)", alpha_hat)),
               each = sum(keep))
)

ggplot(sc_df, aes(truth, fitted)) +
  geom_point(alpha = 0.25, size = 0.6, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ model, nrow = 2) +
  labs(title = "Predicted vs. true albedo",
       x = "True", y = "Posterior mean") +
  theme_minimal(base_size = 12)
