# =============================================================================
# Experiment 1 (SIF): Landcover-aware SAR prior on real OCO-2 SIF data
#
# Same LC-aware W as the albedo experiment, but applied to the real SIF757
# observations. No ground truth exists, so evaluation is:
#   - CV score comparison (standard W vs LC-aware W)
#   - Posterior map and SE inspection
#   - NDVI correlation as external validation
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

# Real SIF observations (no noise added -- this is the actual measurement)
y_sif <- goebel2026::soundings_augmented$SIF_757nm

# Covariates
water_sounding             <- goebel2026::soundings_augmented$proportion_water
X_fixed                    <- cbind(intercept = 1, water = water_sounding)
q                          <- ncol(X_fixed)

water_grid_p               <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0
X_grid                     <- cbind(intercept = 1, water = water_grid_p)

# NDVI on target grid (external validation only -- not used in fitting)
ndvi_grid <- goebel2026::target_grid$mean_ndvi

# -----------------------------------------------------------------------------
# Landcover-aware W
# -----------------------------------------------------------------------------
W_std <- goebel2026::make_W_matrix(goebel2026::target_grid)

make_W_lc <- function(alpha, water_prop, W_template) {
  W_trip <- Matrix::summary(W_template)

  xi    <- water_prop[W_trip$i]
  xj    <- water_prop[W_trip$j]
  raw_w <- pmax(1 - alpha * abs(xi - xj), 0)

  row_s_tab <- tapply(raw_w, W_trip$i, sum)
  row_s_vec <- numeric(nrow(W_template))
  row_s_vec[as.integer(names(row_s_tab))] <- as.numeric(row_s_tab)

  row_s_safe <- pmax(row_s_vec, .Machine$double.eps)
  norm_w     <- raw_w / row_s_safe[W_trip$i]

  Matrix::sparseMatrix(
    i    = W_trip$i,
    j    = W_trip$j,
    x    = as.numeric(norm_w),
    dims = dim(W_template)
  )
}

make_Q_spatial <- function(alpha) {
  W_lc    <- make_W_lc(alpha, water_grid_p, W_std)
  IminusW <- Matrix::Diagonal(p) - W_lc
  Q_sp    <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  Matrix::drop0(Q_sp)
}

# -----------------------------------------------------------------------------
# Augmented design matrix [A | X_fixed]
# -----------------------------------------------------------------------------
lambda_beta <- 0.01
A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

# -----------------------------------------------------------------------------
# RSR constraint function for CV
#
# A_aug has p spatial columns then q covariate columns.
# The RSR constraint is C = [t(X_train) %*% A_train_spatial | 0_{q x q}],
# enforcing orthogonality between the spatial field and the covariates
# on training observations only -- which is what changes fold to fold.
# -----------------------------------------------------------------------------
rsr_constraint <- function(train_idx, A_train_aug) {
  A_train_spatial <- A_train_aug[, seq_len(p), drop = FALSE]
  C_spatial       <- as.matrix(t(X_fixed[train_idx, , drop = FALSE]) %*% A_train_spatial)
  cbind(C_spatial, matrix(0, q, q))
}

# -----------------------------------------------------------------------------
# Q_fun factories
# -----------------------------------------------------------------------------
Q_aug_fun_lc <- function(theta) {
  alpha     <- theta[["alpha"]]
  Q_sp      <- make_Q_spatial(alpha)
  Q_aug     <- Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

Q_aug_fun_std <- function(theta) {
  Q_sp      <- make_Q_spatial(0.0)
  Q_aug     <- Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))
  log_det_Q <- q * log(lambda_beta)
  list(Q = Q_aug, log_det_Q = log_det_Q)
}

# -----------------------------------------------------------------------------
# Tune standard W (alpha = 0): phi only
# -----------------------------------------------------------------------------
cat("=== Tuning standard W (alpha = 0) ===\n")
tuned_std <- tune_cv(
  y          = y_sif,
  A          = A_aug,
  Q_fun      = Q_aug_fun_std,
  theta_init = numeric(0),
  k          = 10L,
  constraint = rsr_constraint,
  verbose    = TRUE
)
phi_std <- tuned_std$phi
cat(sprintf("Standard W: phi = %.4f   CV-MSE = %.4f\n\n",
            phi_std, tuned_std$value))

# -----------------------------------------------------------------------------
# Tune LC-aware W: alpha and phi jointly
# -----------------------------------------------------------------------------
cat("=== Tuning LC-aware W (alpha in [0,1]) ===\n")
tuned_lc <- tune_cv(
  y          = y_sif,
  A          = A_aug,
  Q_fun      = Q_aug_fun_lc,
  theta_init = c(alpha = 0.5),
  lower      = c(alpha = 0.0),
  upper      = c(alpha = 1.0),
  k          = 10L,
  constraint = rsr_constraint,
  verbose    = TRUE,
  parallel = TRUE
)
alpha_hat <- tuned_lc$theta[["alpha"]]
phi_lc    <- tuned_lc$phi
cat(sprintf("LC-aware W: alpha = %.4f   phi = %.4f   CV-MSE = %.4f\n\n",
            alpha_hat, phi_lc, tuned_lc$value))

cat(sprintf("CV-MSE reduction: %.2f%%\n",
            100 * (tuned_std$value - tuned_lc$value) / tuned_std$value))

# -----------------------------------------------------------------------------
# Fit helper
# -----------------------------------------------------------------------------
make_fit <- function(alpha, phi) {
  Q_sp  <- make_Q_spatial(alpha)
  Q_aug <- Matrix::bdiag(Q_sp, lambda_beta / phi * Matrix::Diagonal(q))
  fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = Q_aug,
    phi    = phi,
    solver = "cholesky"
  )
}

# Full-data RSR constraint (all n observations)
C_aug_full <- cbind(t(X_fixed) %*% A, matrix(0, q, q))

# -----------------------------------------------------------------------------
# Fit: standard W
# -----------------------------------------------------------------------------
fit_std <- make_fit(0.0, phi_std)

beta_hat_std <- fit_std$posterior_mean[p + seq_len(q)]
cat(sprintf("Standard W  beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_std[1], beta_hat_std[2]))

A_pred      <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred_std    <- as.vector(A_pred %*% fit_std$posterior_mean)
se_pred_std <- fastblm::posterior_se(fit_std, A_new = A_pred)

fit_std_rsr      <- fastblm::constrain(fit_std, C_aug_full)
pred_std_rsr     <- as.vector(A_pred %*% fit_std_rsr$posterior_mean)
se_pred_std_rsr  <- fastblm::posterior_se(fit_std_rsr, A_new = A_pred)
beta_hat_std_rsr <- fit_std_rsr$posterior_mean[p + seq_len(q)]
cat(sprintf("Standard W RSR beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_std_rsr[1], beta_hat_std_rsr[2]))

# -----------------------------------------------------------------------------
# Fit: LC-aware W
# -----------------------------------------------------------------------------
fit_lc <- make_fit(alpha_hat, phi_lc)

beta_hat_lc <- fit_lc$posterior_mean[p + seq_len(q)]
cat(sprintf("LC-aware W  beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_lc[1], beta_hat_lc[2]))

pred_lc    <- as.vector(A_pred %*% fit_lc$posterior_mean)
se_pred_lc <- fastblm::posterior_se(fit_lc, A_new = A_pred)

fit_lc_rsr      <- fastblm::constrain(fit_lc, C_aug_full)
pred_lc_rsr     <- as.vector(A_pred %*% fit_lc_rsr$posterior_mean)
se_pred_lc_rsr  <- fastblm::posterior_se(fit_lc_rsr, A_new = A_pred)
beta_hat_lc_rsr <- fit_lc_rsr$posterior_mean[p + seq_len(q)]
cat(sprintf("LC-aware W RSR beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_lc_rsr[1], beta_hat_lc_rsr[2]))

# -----------------------------------------------------------------------------
# NDVI correlation (external validation)
# -----------------------------------------------------------------------------
observed_cells <- colSums(A) > 0
ndvi_valid     <- observed_cells & !is.na(ndvi_grid) & water_grid_p < 0.5

ndvi_cor <- function(label, pred) {
  r2 <- cor(pred[ndvi_valid], ndvi_grid[ndvi_valid])^2
  cat(sprintf("%-40s  NDVI R2 = %.4f  (n=%d cells)\n",
              label, r2, sum(ndvi_valid)))
  r2
}

cat("\n--- NDVI correlation (water < 50%%) ---\n")
ndvi_cor("Standard W — unconstrained",         pred_std)
ndvi_cor("Standard W — RSR",                   pred_std_rsr)
ndvi_cor(sprintf("LC-aware (α=%.3f) — unconstrained", alpha_hat), pred_lc)
ndvi_cor(sprintf("LC-aware (α=%.3f) — RSR",    alpha_hat), pred_lc_rsr)

summary_df <- data.frame(
  Model = c(
    "Standard W",
    "Standard W + RSR",
    sprintf("LC-aware W (α=%.3f)", alpha_hat),
    sprintf("LC-aware W (α=%.3f) + RSR", alpha_hat)
  ),
  CV_MSE  = c(tuned_std$value, NA, tuned_lc$value, NA),
  phi     = c(phi_std, phi_std, phi_lc, phi_lc),
  NDVI_R2 = c(
    ndvi_cor("", pred_std),
    ndvi_cor("", pred_std_rsr),
    ndvi_cor("", pred_lc),
    ndvi_cor("", pred_lc_rsr)
  ),
  beta_water = c(
    beta_hat_std[2], beta_hat_std_rsr[2],
    beta_hat_lc[2],  beta_hat_lc_rsr[2]
  )
)

cat("\n")
knitr::kable(summary_df, digits = 4, row.names = FALSE)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
mask <- function(x) ifelse(observed_cells, x, NA_real_)

plot_sf <- target_proj
plot_sf$ndvi         <- ndvi_grid
plot_sf$water        <- water_grid_p
plot_sf$pred_std     <- mask(pred_std)
plot_sf$pred_std_rsr <- mask(pred_std_rsr)
plot_sf$pred_lc      <- mask(pred_lc)
plot_sf$pred_lc_rsr  <- mask(pred_lc_rsr)
plot_sf$se_std       <- mask(se_pred_std)
plot_sf$se_lc        <- mask(se_pred_lc)
plot_sf$diff_pred    <- mask(pred_lc - pred_std)
plot_sf$diff_se      <- mask(se_pred_lc - se_pred_std)

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

lim_sif  <- range(c(pred_std, pred_lc, pred_lc_rsr), na.rm = TRUE)
lim_se   <- range(c(se_pred_std, se_pred_lc), na.rm = TRUE)
lim_diff <- max(abs(plot_sf$diff_pred), na.rm = TRUE) * c(-1, 1)

(sf_plot(plot_sf, "pred_std",     "Posterior mean — standard W",
         lims = lim_sif, midpt = mean(lim_sif)) |
    sf_plot(plot_sf, "pred_lc",      sprintf("Posterior mean — LC-aware W (α=%.2f)", alpha_hat),
            lims = lim_sif, midpt = mean(lim_sif)) |
    sf_plot(plot_sf, "pred_lc_rsr",  sprintf("Posterior mean — LC-aware RSR (α=%.2f)", alpha_hat),
            lims = lim_sif, midpt = mean(lim_sif))) +
  plot_annotation(
    title = sprintf("Downscaled SIF757   phi_std=%.3f  phi_lc=%.3f  alpha=%.3f",
                    phi_std, phi_lc, alpha_hat)
  )

sf_plot(plot_sf, "diff_pred",
        sprintf("Δ posterior mean (LC-aware − standard, α=%.2f)", alpha_hat),
        lims = lim_diff, midpt = 0)

(sf_plot(plot_sf, "se_std", "SE — standard W",
         lo = "#fff7ec", mi = "#fc8d59", hi = "#7f0000",
         midpt = mean(lim_se), lims = lim_se) |
    sf_plot(plot_sf, "se_lc",  sprintf("SE — LC-aware W (α=%.2f)", alpha_hat),
            lo = "#fff7ec", mi = "#fc8d59", hi = "#7f0000",
            midpt = mean(lim_se), lims = lim_se))

(sf_plot(plot_sf, "ndvi", "NDVI",
         lo = "#d73027", mi = "#ffffbf", hi = "#1a9850", midpt = 0.5) |
    sf_plot(plot_sf, "pred_lc_rsr",
            sprintf("SIF — LC-aware RSR (α=%.2f)", alpha_hat),
            lims = lim_sif, midpt = mean(lim_sif)))

keep_ndvi <- ndvi_valid
sc_df <- data.frame(
  ndvi   = rep(ndvi_grid[keep_ndvi], 2),
  sif    = c(pred_std_rsr[keep_ndvi], pred_lc_rsr[keep_ndvi]),
  model  = rep(c("Standard W + RSR",
                 sprintf("LC-aware RSR (α=%.2f)", alpha_hat)),
               each = sum(keep_ndvi)),
  water  = rep(water_grid_p[keep_ndvi], 2)
)

ggplot(sc_df, aes(ndvi, sif)) +
  geom_point(aes(colour = water), alpha = 0.4, size = 0.8) +
  scale_colour_gradient(low = "#2c7bb6", high = "#d7191c",
                        name = "Water\nfraction") +
  geom_smooth(method = "lm", se = FALSE, colour = "black",
              linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ model) +
  labs(title = "Downscaled SIF vs NDVI (water < 50%)",
       x = "NDVI", y = "Posterior mean SIF757") +
  theme_minimal(base_size = 12)
