## =============================================================================
## example_diagnostics.R
##
## Run this after ablations_1_3.R — assumes all fit objects are in the
## workspace. Source diagnostics.R first.
## =============================================================================

#source("diagnostics.R")   # or wherever you saved it
library(patchwork)

# =============================================================================
# Step 1: Build a diagnostics object for each model
# =============================================================================
# blm_diagnostics(fit, se, truth, grid_sf, label, C_check = NULL)
# grid_sf must be an sf with p rows (one per latent grid cell)

diag_baseline <- blm_diagnostics(
  fit     = fit_cv,
  se      = se_cv,
  truth   = truth,
  grid_sf = target_proj,          # sf object, p rows
  label   = "Baseline (random CV)"
)

diag_spatial <- blm_diagnostics(
  fit     = fit_cv_spatial,
  se      = se_cv_spatial,
  truth   = truth,
  grid_sf = target_proj,
  label   = "Spatial k-means CV"
)

diag_rsr <- blm_diagnostics(
  fit     = fit_rsr,
  se      = se_rsr,
  truth   = truth,
  grid_sf = target_proj,
  label   = "RSR (constrained + covariates)",
  C_check = C_global              # verifies ||C mu||_inf ~ 0
)

# Print scalar summaries for any individual model
print(diag_baseline)
print(diag_rsr)

# =============================================================================
# Step 2: Comparison table across all models
# =============================================================================

comp <- compare_diagnostics(list(
  baseline = diag_baseline,
  spatial  = diag_spatial,
  rsr      = diag_rsr
))
print(comp, digits = 4)

# =============================================================================
# Step 3: Four-panel spatial map for a single model
# =============================================================================
# Posterior mean | Posterior SE
# Error          | 95% CI coverage

plot(diag_rsr)

# Save to file if needed
ggsave("rsr_spatial_maps.png", plot = plot(diag_rsr),
       width = 10, height = 7, dpi = 150)

# =============================================================================
# Step 4: Scatter truth vs posterior mean (like paper Fig 4)
# =============================================================================
# Colour by water fraction to check covariate separation

water_grid <- goebel2026::target_grid$prop_water
water_grid[is.na(water_grid)] <- 0

scatter_truth_vs_mean(diag_baseline, colour_var = water_grid,
                      colour_label = "Prop. water")

scatter_truth_vs_mean(diag_rsr, colour_var = water_grid,
                      colour_label = "Prop. water")

# Side-by-side comparison of baseline vs RSR scatter
p_scatter <- scatter_truth_vs_mean(diag_baseline, colour_var = water_grid,
                                   colour_label = "Prop. water") |
  scatter_truth_vs_mean(diag_rsr,      colour_var = water_grid,
                        colour_label = "Prop. water")
p_scatter

# =============================================================================
# Step 5: CV optimisation trace  (one per tuned object)
# =============================================================================
# Shows rho evaluations as points coloured by profiled log(phi),
# red x marks the optimum

plot_cv_profile(tuned_cv_spatial)   # spatial folds run
plot_cv_profile(tuned_cv_rsr)       # RSR run

# =============================================================================
# Step 6: Fold assignment + per-fold errors on sounding map
# =============================================================================
# Three panels: fold assignments | per-fold RMSE | pointwise sounding residuals

plot_fold_errors(
  soundings_sf = soundings_proj,
  folds        = spatial_folds,
  y            = y_alb,
  A            = A,
  fit          = fit_cv_spatial
)

# Compare random vs spatial fold assignments side by side
# (make random folds first if you haven't already)
set.seed(42)
random_folds <- sample(rep_len(seq_len(10L), length(y_alb)))

p_folds_random  <- plot_fold_errors(soundings_proj, random_folds,
                                    y_alb, A, fit_cv)
p_folds_spatial <- plot_fold_errors(soundings_proj, spatial_folds,
                                    y_alb, A, fit_cv_spatial)

# Print separately — patchwork nesting gets unwieldy for 3+3 panels
p_folds_random
p_folds_spatial

# =============================================================================
# Step 7: Cross-model map comparisons on a shared colour scale
# =============================================================================

all_diags <- list(
  "Baseline"   = diag_baseline,
  "Spatial CV" = diag_spatial,
  "RSR"        = diag_rsr
)

# Posterior means side by side
plot_compare_maps(all_diags, field = "posterior_mean")

# Errors side by side — same diverging scale makes differences obvious
plot_compare_maps(all_diags, field = "error")

# Posterior SE — shows where RSR changes uncertainty vs baseline
plot_compare_maps(all_diags, field = "posterior_se")

# CI coverage map — look for spatially structured miscalibration
plot_compare_maps(all_diags, field = "in_ci")

# =============================================================================
# Step 8: RSR constraint check
# =============================================================================
# ||C * posterior_mean||_inf should be near machine epsilon for a
# correctly constrained fit. Printed by blm_diagnostics but check directly:

cat(sprintf("RSR constraint residual ||C mu||_inf = %.2e\n",
            max(abs(as.numeric(C_global %*% fit_rsr$posterior_mean)))))

# For comparison, the unconstrained fit:
cat(sprintf("Unconstrained           ||C mu||_inf = %.2e\n",
            max(abs(as.numeric(C_global %*% fit_rsr_unconstrained$posterior_mean)))))
