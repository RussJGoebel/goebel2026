# diagnostic_water_zero.R
#
# Self-contained diagnostic: fit only on majority-land pixels,
# drop soundings that only overlap majority-water pixels,
# then zero out water pixels in the final result.

library(fastblm)
library(goebel2026)
library(Matrix)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

d_shared           <- goebel2026::setup_shared
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_latent_water     <- d_shared$X_latent_water
fine_grid_buffered <- d_shared$fine_grid_buffered

y_sif      <- goebel2026::soundings_augmented$SIF_757nm
target_idx <- which(fine_grid_buffered$n_intersects > 0)
p          <- ncol(A)

land_frac  <- 1 - as.numeric(X_latent_water[, 2])
land_mask  <- land_frac >= 0.5
water_mask <- !land_mask

message(sprintf("Total pixels: %d  |  majority land: %d  |  majority water: %d",
                p, sum(land_mask), sum(water_mask)))

# ------------------------------------------------------------------------------
# 2. Subset A columns to land pixels, drop all-water soundings
# ------------------------------------------------------------------------------

A_land <- A[, land_mask, drop = FALSE]

# Drop soundings with zero overlap on any land pixel
sounding_has_land <- Matrix::rowSums(abs(A_land)) > 0
A_filt            <- A_land[sounding_has_land, , drop = FALSE]
y_filt            <- y_sif[sounding_has_land]

message(sprintf("Soundings: total=%d  kept=%d  dropped (all-water)=%d",
                length(y_sif), sum(sounding_has_land), sum(!sounding_has_land)))

# ------------------------------------------------------------------------------
# 3. Build SAR prior on land pixels only
#
# Re-normalize W so rows still sum correctly after dropping water neighbors.
# ------------------------------------------------------------------------------

p_land      <- ncol(A_filt)
W_land      <- W_queen[land_mask, land_mask]
row_sums    <- Matrix::rowSums(W_land)
row_sums[row_sums == 0] <- 1  # isolated pixels
W_land_norm <- Matrix::Diagonal(x = 1 / row_sums) %*% W_land
S_land      <- Matrix::Diagonal(p_land) - W_land_norm
Q_land      <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_land)))
# Small ridge for numerical stability -- some boundary land pixels have
# few land neighbors after dropping water adjacencies, making Q near-singular
Q_land      <- Q_land + Matrix::Diagonal(p_land) * 1e-6

message(sprintf("Land-only system: p=%d  n_snd=%d  nnz(Q)=%d",
                p_land, nrow(A_filt), Matrix::nnzero(Q_land)))

# ------------------------------------------------------------------------------
# 4. Fit
# ------------------------------------------------------------------------------

PHI <- 5
message(sprintf("\nFitting at phi=%d...", PHI))

fit_land <- fastblm::fit_fastblm(
  y      = y_filt,
  A      = A_filt,
  Q      = Q_land,
  phi    = 1/PHI,
  solver = "cholesky"
)

message(sprintf("  sigma2e=%.4g", fit_land$sigma2e))

# ------------------------------------------------------------------------------
# 5. Fill full grid: land = fitted, water = 0
# ------------------------------------------------------------------------------

mu_full             <- numeric(p)
mu_full[land_mask]  <- fit_land$posterior_mean
mu_full[water_mask] <- 0

message(sprintf("\nPosterior mean at land pixels:  mean=%.4f  max=%.4f",
                mean(mu_full[land_mask]), max(mu_full[land_mask])))
message(sprintf("Posterior mean at water pixels: all exactly 0 (n=%d)", sum(water_mask)))

# ------------------------------------------------------------------------------
# 6. Scatter and map
# ------------------------------------------------------------------------------

grid_diag            <- goebel2026::target_grid
grid_diag            <- grid_diag[grid_diag$n_intersects > 0, ]
grid_diag$mu         <- mu_full[target_idx]
grid_diag$water_frac <- 1 - land_frac[target_idx]

p_scatter <- ggplot(grid_diag, aes(x = water_frac, y = mu)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red") +
  labs(x = "Water fraction", y = "Posterior mean SIF",
       title = sprintf("Water-zero: land-only fit + zeroed water pixels, phi=%d", PHI),
       subtitle = sprintf("majority-water in target: n=%d  all zero  |  land mean=%.3f",
                          sum(grid_diag$water_frac >= 0.5),
                          mean(grid_diag$mu[grid_diag$water_frac < 0.5]))) +
  theme_minimal()

print(p_scatter)

p_map <- ggplot(grid_diag) +
  geom_sf(aes(fill = mu*(mu <= 3)), color = NA) +
  scale_fill_viridis_c(name = "SIF") +
  labs(title = sprintf("Water-zero: land-only fit + zeroed water, phi=%d", PHI)) +
  theme_minimal()

print(p_map)

message("\nDone.")
