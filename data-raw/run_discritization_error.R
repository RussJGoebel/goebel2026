# data-raw/run_discretization_error.R
#
# Estimates the relative operator error of the A matrix discretization at
# multiple latent grid resolutions, compared to the ground truth integration
# provided by exactextractr (soundings_augmented$mean_albedo).
#
# The quantity estimated is the relative error on the albedo field:
#
#   rel_error = ||A_approx * y_latent - A_true * y_latent|| / ||A_true * y_latent||
#
# where:
#   A_true * y_latent  = soundings_augmented$mean_albedo  (exactextractr ground truth)
#   A_approx * y_latent = A %*% y_latent_at_resolution    (sparse matrix approximation)
#
# This is a lower bound on the operator norm of the error matrix E = A_approx - A_true,
# evaluated on a single test vector (the albedo field). Results are saved as
# results_discretization_error.
#
# Resolutions tested: 330m (current), 165m, 110m
# Each requires building a new target grid and A matrix -- no model fitting.
#
# Output saved via usethis::use_data():
#   results_discretization_error  -- data frame with resolution, n_cells,
#                                    rel_error, rmse, r2 per resolution

library(goebel2026)
library(spatintegrate)
library(terra)
library(sf)
library(future)
library(future.apply)
library(usethis)

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load shared objects
# ------------------------------------------------------------------------------

message("Loading shared objects...")

keep_idx       <- goebel2026::setup_sif$keep_idx
soundings_proj <- spatintegrate::ensure_projected(
  goebel2026::soundings_augmented
)[keep_idx, ]

# Ground truth: direct 10m->coarse integration via exactextractr
y_true <- goebel2026::soundings_augmented$mean_albedo[keep_idx]
message(sprintf("  n soundings: %d  y_true mean=%.4f  sd=%.4f",
                length(y_true), mean(y_true, na.rm = TRUE), sd(y_true, na.rm = TRUE)))

# Raster for albedo extraction at each resolution
median_albedo_rast <- terra::rast(
  "/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif"
)

# Denominator for relative error (fixed across resolutions)
norm_y_true <- sqrt(sum(y_true^2, na.rm = TRUE))

# ------------------------------------------------------------------------------
# 2. Helper: build grid, aggregate albedo, build A, compute errors
# ------------------------------------------------------------------------------

run_resolution <- function(res_m) {
  message(sprintf("\n== Resolution: %dm ==", res_m))

  buf <- 25 * res_m

  # Build target grid at this resolution
  message("  Building target grid...")
  tg <- spatintegrate::make_square_grid_in_crs(
    soundings_proj,
    res_m,
    buffer = buf
  )

  # Reproject to raster CRS for extraction
  tg_rast_crs <- sf::st_transform(tg, sf::st_crs(median_albedo_rast))

  # Aggregate 10m albedo to this resolution
  message("  Aggregating albedo raster...")
  tg_rast_crs <- spatintegrate::summarize_raster_mean_over_polygons(
    median_albedo_rast,
    tg_rast_crs,
    stats_col_prefix = "mean_"
  )

  # Handle column name (summarize_raster appends raster layer name)
  alb_col <- grep("mean_", names(tg_rast_crs), value = TRUE)[1]
  y_latent <- sf::st_drop_geometry(tg_rast_crs)[[alb_col]]
  message(sprintf("  n_cells=%d  y_latent mean=%.4f  NAs=%d",
                  length(y_latent), mean(y_latent, na.rm = TRUE), sum(is.na(y_latent))))

  # Replace NAs with 0 (edge cells outside raster extent)
  y_latent[is.na(y_latent)] <- 0

  # Reproject grid back to soundings CRS for A construction
  tg_snd_crs <- sf::st_transform(tg_rast_crs, sf::st_crs(soundings_proj))

  # Build A matrix
  message("  Building A matrix...")
  t_start <- proc.time()
  A <- spatintegrate::compute_overlap_fractions(
    soundings = soundings_proj,
    fine_grid = tg_snd_crs,
    parallel  = TRUE
  )
  t_elapsed <- (proc.time() - t_start)[["elapsed"]]
  message(sprintf("  A: %d x %d  nnz=%d  time=%.1fs",
                  nrow(A), ncol(A), Matrix::nnzero(A), t_elapsed))

  # Approximate integration
  y_approx <- as.numeric(A %*% y_latent)

  # Errors
  err      <- y_approx - y_true
  rel_err  <- sqrt(sum(err^2, na.rm = TRUE)) / norm_y_true
  rmse     <- sqrt(mean(err^2, na.rm = TRUE))
  r2       <- cor(y_approx, y_true, use = "complete.obs")^2

  message(sprintf("  rel_error=%.6f  RMSE=%.6f  R2=%.6f", rel_err, rmse, r2))

  list(
    resolution = res_m,
    n_cells    = nrow(tg),
    nnz_A      = Matrix::nnzero(A),
    rel_error  = rel_err,
    rmse       = rmse,
    r2         = r2,
    build_time = t_elapsed
  )
}

# ------------------------------------------------------------------------------
# 3. Run across resolutions
# ------------------------------------------------------------------------------

resolutions <- c(330, 165, 110)

results_list <- lapply(resolutions, run_resolution)

results_discretization_error <- do.call(rbind, lapply(results_list, as.data.frame))
rownames(results_discretization_error) <- NULL

message("\n== Summary ==")
print(results_discretization_error, digits = 6)

# ------------------------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------------------------

usethis::use_data(results_discretization_error, overwrite = TRUE)

future::plan(future::sequential)
message("\nrun_discretization_error.R complete.")
message("  results_discretization_error saved to data/")
