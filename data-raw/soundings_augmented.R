## =============================================================================
## compute_soundings_augmented.R
##
## Builds soundings_augmented from goebel2026::soundings by attaching:
##   - proportion_water  (from ESA WorldCover, class 6 binary mask)
##   - mean_ndvi         (from Sentinel-2 median NDVI, June 2022)
##   - mean_albedo       (from Sentinel-2 median albedo, June 2022)
##
## Uses exactextractr::exact_extract for all three -- handles overlapping
## sounding footprints correctly (each polygon processed independently).
## Rasters are cropped to soundings extent before extraction to manage memory.
##
## Requires: module load geos before starting R
## =============================================================================

library(goebel2026)
library(spatintegrate)
library(terra)
library(sf)
library(exactextractr)

# =============================================================================
# 1. Load rasters
# =============================================================================

landcover_rast  <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif")
ndvi_rast       <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_ndvi_june_2022.tif")
albedo_rast     <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif")

# =============================================================================
# 2. Project soundings to each raster's CRS for extraction, then crop rasters
# =============================================================================

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)

# Helper: reproject soundings to a raster's CRS and return its extent
soundings_in_crs <- function(rast) {
  sf::st_transform(soundings_proj, sf::st_crs(rast))
}

cat("Cropping rasters to soundings extent...\n")
soundings_lc     <- soundings_in_crs(landcover_rast)
soundings_ndvi   <- soundings_in_crs(ndvi_rast)
soundings_albedo <- soundings_in_crs(albedo_rast)

landcover_crop  <- terra::crop(landcover_rast, terra::ext(terra::vect(soundings_lc)))
ndvi_crop       <- terra::crop(ndvi_rast,      terra::ext(terra::vect(soundings_ndvi)))
albedo_crop     <- terra::crop(albedo_rast,    terra::ext(terra::vect(soundings_albedo)))

cat(sprintf("Landcover crop:  %d x %d\n", nrow(landcover_crop),  ncol(landcover_crop)))
cat(sprintf("NDVI crop:       %d x %d\n", nrow(ndvi_crop),       ncol(ndvi_crop)))
cat(sprintf("Albedo crop:     %d x %d\n", nrow(albedo_crop),     ncol(albedo_crop)))

# =============================================================================
# 3. Build binary water mask (class 6 = permanent water bodies)
# =============================================================================

cat("Building binary water mask...\n")
water_rast <- terra::ifel(landcover_crop == 6, 1, 0)

# =============================================================================
# 4. Extract covariates per sounding via exactextractr
# =============================================================================

cat("Extracting proportion_water...\n")
system.time({
  proportion_water <- exactextractr::exact_extract(water_rast,  soundings_lc,     fun = "mean")
})

cat("Extracting mean_ndvi...\n")
system.time({
  mean_ndvi        <- exactextractr::exact_extract(ndvi_crop,   soundings_ndvi,   fun = "mean")
})

cat("Extracting mean_albedo...\n")
system.time({
  mean_albedo      <- exactextractr::exact_extract(albedo_crop, soundings_albedo, fun = "mean")
})

# =============================================================================
# 5. Clean up: clamp water to [0,1], replace NAs with 0 for water only
# =============================================================================

proportion_water <- pmin(pmax(replace(proportion_water, is.na(proportion_water), 0), 0), 1)

cat("\n--- Summaries ---\n")
cat(sprintf("proportion_water: mean=%.4f  sd=%.4f  prop_zero=%.3f\n",
            mean(proportion_water), sd(proportion_water), mean(proportion_water == 0)))
cat(sprintf("mean_ndvi:        mean=%.4f  sd=%.4f  NAs=%d\n",
            mean(mean_ndvi, na.rm=TRUE), sd(mean_ndvi, na.rm=TRUE), sum(is.na(mean_ndvi))))
cat(sprintf("mean_albedo:      mean=%.4f  sd=%.4f  NAs=%d\n",
            mean(mean_albedo, na.rm=TRUE), sd(mean_albedo, na.rm=TRUE), sum(is.na(mean_albedo))))

# =============================================================================
# 6. Attach to soundings and save as soundings_augmented
# =============================================================================

soundings_augmented <- goebel2026::soundings
soundings_augmented$proportion_water <- proportion_water
soundings_augmented$mean_ndvi        <- mean_ndvi
soundings_augmented$mean_albedo      <- mean_albedo

usethis::use_data(soundings_augmented, overwrite = TRUE)
cat("Saved soundings_augmented.\n")
