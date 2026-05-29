# data-raw/diagnostic_spatial_residuals.R
#
# Spatial residual diagnostic for the canonical SIF fit.
#
# Depends on setup.R having already been run (outlier removed there).
# All data objects loaded from the goebel2026 package.
#
# PURPOSE
# -------
# Assesses whether there is residual spatial autocorrelation in the fitted
# model after accounting for the observation geometry. A correctly-specified
# model should produce standardized prediction residuals that are iid N(0,1)
# with no spatial structure.
#
# KEY IDEA
# --------
# Raw residuals are not informative about spatial misspecification because
# sounding informativeness varies strongly with observation density. We
# therefore use standardized prediction residuals:
#
#   r_i = (y_i - A_i mu) / sqrt(Var_pred_i)
#
# where Var_pred_i = sigma2e * (A_i K^{-1} A_i') + sigma2e / [R^{-1}]_ii
#
#   - sigma2e * A_i K^{-1} A_i'  is the posterior variance of the fitted value
#   - sigma2e / [R^{-1}]_ii      is the per-sounding observation noise variance
#
# Spatial structure in the variogram of r_i is evidence of misspecification.
#
# OUTPUT
# ------
#   diagnostic_spatial_residuals.pdf  -- four-panel figure:
#     A: spatial map of standardized residuals
#     B: std. residual vs. predictive SD (information content check)
#     C: empirical variogram of std. residuals (main diagnostic)
#     D: Q-Q plot

# ==============================================================================
# 0. Packages
# ==============================================================================

required <- c("fastblm", "goebel2026", "Matrix",
              "gstat", "sf", "ggplot2", "patchwork")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0L)
  stop("Please install: ", paste(missing, collapse = ", "))

library(fastblm)
library(goebel2026)
library(Matrix)
library(gstat)
library(sf)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. Load data (outlier already removed in setup.R)
# ==============================================================================

message("Loading data...")

res      <- goebel2026::results_sif_canonical
d_shared <- goebel2026::setup_shared
d_sif    <- goebel2026::setup_sif

W_queen     <- d_shared$W_queen
X_obs_water <- d_shared$X_obs_water  # consistent with y_sif -- filter in section 1 of setup.R

y_sif     <- d_sif$y
R_inv_sif <- d_sif$R_inv
A         <- d_shared$A_flat         # built after filter, so already n-1 rows

sounding_coords <- as.data.frame(
  sf::st_coordinates(sf::st_centroid(d_shared$soundings_proj))
)
names(sounding_coords) <- c("x", "y")

p           <- ncol(A)
q           <- ncol(X_obs_water)
lambda_beta <- 0.01
n           <- length(y_sif)

message(sprintf("n soundings = %d, p latent cells = %d", n, p))

rho_hat <- res$rho_opt
phi_hat <- res$phi
message(sprintf("Tuned hyperparameters from canonical fit: rho=%.4f  phi=%.4f",
                rho_hat, phi_hat))

# ==============================================================================
# 2. Rebuild model objects at tuned hyperparameters
# ==============================================================================

S     <- Matrix::Diagonal(nrow(W_queen)) - rho_hat * W_queen
Q_sp  <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
Q_aug <- Matrix::drop0(Matrix::forceSymmetric(
  Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))
))

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

C_aug <- cbind(
  as.matrix(t(X_obs_water) %*% A),
  matrix(0, nrow = q, ncol = q)
)

# ==============================================================================
# 3. Fit at tuned hyperparameters
# ==============================================================================

message("Fitting at tuned hyperparameters...")
fit <- fastblm::fit_fastblm(
  y      = y_sif,
  A      = A_aug,
  Q      = Q_aug,
  phi    = phi_hat,
  R_inv  = R_inv_sif,
  solver = "cholesky"
)
fit <- fastblm::constrain(fit, C_aug)
message(sprintf("  sigma2e = %.4g", fit$sigma2e))

# ==============================================================================
# 4. Standardized prediction residuals
# ==============================================================================

message("Computing standardized residuals...")

fitted_vals    <- as.numeric(A_aug %*% fit$posterior_mean)
raw_resid      <- y_sif - fitted_vals

# Posterior SD of fitted value A_i mu
post_sd_fitted <- fastblm::posterior_se(fit, A_new = A_aug)

# Per-sounding observation noise SD from the diagonal of R^{-1}
R_inv_diag     <- Matrix::diag(R_inv_sif)
obs_noise_var  <- fit$sigma2e / R_inv_diag

# Total predictive variance = posterior variance of fitted + obs noise
pred_var  <- post_sd_fitted^2 + obs_noise_var
pred_sd   <- sqrt(pred_var)
std_resid <- raw_resid / pred_sd

message(sprintf(
  "Standardized residuals: mean = %.3f (expect ~0),  sd = %.3f (expect ~1)",
  mean(std_resid), sd(std_resid)
))
message(sprintf("  |r| > 2: %.1f%%  (expect ~5%%)",
                100 * mean(abs(std_resid) > 2)))
message(sprintf("  |r| > 3: %.1f%%  (expect ~0.3%%)",
                100 * mean(abs(std_resid) > 3)))

# ==============================================================================
# 5. Empirical variogram of standardized residuals
# ==============================================================================

message("Computing empirical variogram of standardized residuals...")

df <- data.frame(
  x         = sounding_coords$x,
  y         = sounding_coords$y,
  std_resid = std_resid,
  pred_sd   = pred_sd
)

pts_sf <- sf::st_as_sf(df, coords = c("x", "y"), crs = 32619L)
sp_df  <- as(pts_sf, "Spatial")

bbox     <- sf::st_bbox(pts_sf)
max_dist <- 0.5 * sqrt((bbox["xmax"] - bbox["xmin"])^2 +
                         (bbox["ymax"] - bbox["ymin"])^2)

vgm_emp <- gstat::variogram(
  std_resid ~ 1,
  data   = sp_df,
  cutoff = max_dist,
  width  = max_dist / 15L
)

vgm_fit <- tryCatch(
  gstat::fit.variogram(
    vgm_emp,
    gstat::vgm(
      psill  = var(std_resid),
      model  = "Exp",
      range  = max_dist / 3,
      nugget = 0
    )
  ),
  error = function(e) {
    message("  Note: variogram model fit failed, showing empirical only")
    NULL
  }
)

message("Variogram (first 5 lags):")
print(head(vgm_emp[, c("dist", "gamma", "np")], 5))
message("  Expected sill under no autocorrelation: 1.0")

if (!is.null(vgm_fit)) {
  message("Fitted variogram model:")
  print(vgm_fit)
}

# ==============================================================================
# 6. Four-panel figure
# ==============================================================================

message("Generating figure...")

clamp <- function(x, lo = -3.5, hi = 3.5) pmax(pmin(x, hi), lo)

p_map <- ggplot(df, aes(x = x, y = y, colour = clamp(std_resid))) +
  geom_point(size = 1.0, alpha = 0.75) +
  scale_colour_distiller(
    palette = "RdBu",
    limits  = c(-3.5, 3.5),
    name    = "Std.\nresidual\n(clamped\n±3.5)"
  ) +
  coord_equal() +
  labs(
    title    = "A: Standardized prediction residuals",
    subtitle = "No spatial clustering expected under correct specification",
    x = "Easting (m)", y = "Northing (m)"
  ) +
  theme_bw(base_size = 10)

p_sd <- ggplot(df, aes(x = pred_sd, y = std_resid)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_hline(
    yintercept = c(-2, 0, 2),
    linetype   = c("dashed", "solid", "dashed"),
    colour     = "steelblue"
  ) +
  geom_smooth(method = "loess", se = TRUE,
              colour = "tomato", linewidth = 0.8) +
  labs(
    title    = "B: Std. residual vs. predictive SD",
    subtitle = "Loess should be flat; slope indicates heteroskedastic misspecification",
    x = "Posterior predictive SD", y = "Standardized residual"
  ) +
  theme_bw(base_size = 10)

p_vgm <- ggplot(vgm_emp, aes(x = dist, y = gamma)) +
  geom_point(aes(size = np), colour = "grey30") +
  scale_size_continuous(name = "Pairs", range = c(1, 4)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "steelblue", linewidth = 0.8) +
  labs(
    title    = "C: Variogram of standardized residuals",
    subtitle = "Flat at sill = 1 indicates no residual spatial autocorrelation",
    x = "Distance (m)", y = "Semivariance"
  ) +
  theme_bw(base_size = 10)

if (!is.null(vgm_fit)) {
  vgm_line <- gstat::variogramLine(vgm_fit, maxdist = max_dist, n = 200L)
  p_vgm <- p_vgm +
    geom_line(data = vgm_line, aes(x = dist, y = gamma),
              colour = "tomato", linewidth = 0.8, inherit.aes = FALSE)
}

p_qq <- ggplot(df, aes(sample = std_resid)) +
  stat_qq(alpha = 0.35, size = 0.7) +
  stat_qq_line(colour = "tomato", linewidth = 0.8) +
  labs(
    title    = "D: Q-Q plot",
    subtitle = "Should follow N(0,1) under correct specification",
    x = "Theoretical quantiles", y = "Sample quantiles"
  ) +
  theme_bw(base_size = 10)

combined <- (p_map | p_sd) / (p_vgm | p_qq) +
  patchwork::plot_annotation(
    title    = "Spatial residual diagnostic — canonical SIF fit",
    subtitle = sprintf(
      paste0("n = %d soundings  |  rho = %.3f  |  phi = %.3f  |  ",
             "mean(r) = %.3f  |  sd(r) = %.3f  |  ",
             "variogram nugget = %.3f  |  sill = %.3f  |  range = %.0f m"),
      n, rho_hat, phi_hat,
      mean(std_resid), sd(std_resid),
      if (!is.null(vgm_fit)) vgm_fit$psill[1] else NA,
      if (!is.null(vgm_fit)) sum(vgm_fit$psill)  else NA,
      if (!is.null(vgm_fit)) vgm_fit$range[2]    else NA
    )
  )

ggsave("diagnostic_spatial_residuals.pdf", combined,
       width = 12, height = 9, device = "pdf")
message("Saved: diagnostic_spatial_residuals.pdf")

# ==============================================================================
# 7. Save results object to package for use in supplement
# ==============================================================================

message("Saving diagnostic results to package...")

results_diagnostic_residuals <- list(
  n_soundings  = n,
  rho          = rho_hat,
  phi          = phi_hat,
  sigma2e      = fit$sigma2e,
  std_resid    = std_resid,
  pred_sd      = pred_sd,
  raw_resid    = raw_resid,
  vgm_emp      = vgm_emp,
  vgm_fit      = vgm_fit,
  resid_mean   = mean(std_resid),
  resid_sd     = sd(std_resid),
  pct_gt2      = mean(abs(std_resid) > 2),
  pct_gt3      = mean(abs(std_resid) > 3),
  coords       = sounding_coords
)

usethis::use_data(results_diagnostic_residuals, overwrite = TRUE)
message("Saved: results_diagnostic_residuals")

message("\n=== Final summary ===")
message(sprintf("  n soundings:        %d", n))
message(sprintf("  Residual mean:      %.4f  (expect ~0)", mean(std_resid)))
message(sprintf("  Residual sd:        %.4f  (expect ~1)", sd(std_resid)))
message(sprintf("  |r| > 2:            %.1f%%  (expect ~5%%)",
                100 * mean(abs(std_resid) > 2)))
if (!is.null(vgm_fit)) {
  message(sprintf("  Variogram nugget:   %.4f", vgm_fit$psill[1]))
  message(sprintf("  Variogram sill:     %.4f  (expect ~1)", sum(vgm_fit$psill)))
  message(sprintf("  Variogram range:    %.0f m", vgm_fit$range[2]))
}
