## =============================================================================
## reml_diagnostics.R
##
## Runs REML on the synthetic albedo data (mirroring the baseline example),
## then produces the full diagnostic suite.
##
## Dependencies: fastblm, spatintegrate, goebel2026, Matrix, sf,
##               ggplot2, patchwork
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)

source("diagnostics.R")

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

W <- goebel2026::make_W_matrix(goebel2026::target_grid)

Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  list(Q = Matrix::drop0(Q))
}

truth <- goebel2026::target_grid$mean_albedo

# =============================================================================
# 2. REML tuning
# =============================================================================

tuned_reml <- fastblm::tune_reml(
  y             = y_alb,
  A             = A,
  Q_fun         = Q_fun,
  theta_init    = c(rho = 0.9),
  lower         = 0.01,
  upper         = 0.999,
  solver        = "cholesky",
  logdet_method = "cholesky",
  verbose       = TRUE
)

cat(sprintf("\nREML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4f\n",
            tuned_reml$theta[["rho"]], tuned_reml$phi, tuned_reml$sigma2e))

# =============================================================================
# 3. Fit at tuned hyperparameters
# =============================================================================

rho_opt <- tuned_reml$theta[["rho"]]
Q_opt   <- Matrix::forceSymmetric(
  Matrix::crossprod(Matrix::Diagonal(nrow(W)) - rho_opt * W)
)

fit_reml <- fastblm::fit_fastblm(
  y_alb, A, Q_opt, tuned_reml$phi, solver = "cholesky"
)

se_reml <- fastblm::posterior_se(fit_reml)

# =============================================================================
# 4. Save ablation and run diagnostics
# =============================================================================

water_grid <- goebel2026::target_grid$prop_water
water_grid[is.na(water_grid)] <- 0

save_ablation(
  model    = list(
    fit          = fit_reml,
    se           = se_reml,
    tuned        = tuned_reml,
    colour_var   = water_grid,
    colour_label = "Prop. water"
  ),
  run_name = "reml_nocov",
  tags     = list(
    tuning      = "reml",
    covariates  = "none",
    constraint  = "none",
    folds       = NA
  ),
  truth   = truth,
  grid_sf = target_proj,
  A       = A
)

# Later, to reload and compare with other runs:
#
#   diag_list <- load_ablations(filter = list(tuning = "reml"))
#   compare_ablations(diag_list, tag = "reml_variants")
#
#   # or load specific runs by name:
#   diag_list <- load_ablations(run_names = c("reml_nocov", "spatial_cv_water_rsr"))
#   compare_ablations(diag_list, tag = "reml_vs_spatial_rsr")
