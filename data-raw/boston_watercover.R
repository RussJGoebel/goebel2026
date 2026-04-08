## =============================================================================
## compute_sounding_water_v7.R
##
## Computes proportion_water for soundings via exactextractr::exact_extract.
## Handles overlapping soundings correctly and is fast on large rasters.
##
## Requires: module load geos before starting R
## =============================================================================

library(goebel2026)
library(spatintegrate)
library(terra)
library(sf)
library(exactextractr)

# =============================================================================
# 1. Load data
# =============================================================================

landcover_rast <- terra::rast(
  "/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif"
)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)
soundings_wgs  <- sf::st_transform(soundings_proj, sf::st_crs(landcover_rast))

# =============================================================================
# 2. Binary water mask (0/1, no NAs -- mean = fraction water)
# =============================================================================

cat("Building binary water mask (class 6)...\n")
is_water <- terra::ifel(landcover_rast == 6, 1, 0)

cat(sprintf("Raster size: %d x %d = %.1fM pixels\n",
            nrow(is_water), ncol(is_water),
            nrow(is_water) * ncol(is_water) / 1e6))

# =============================================================================
# 3. Extract water fraction per sounding via exactextractr
# =============================================================================

cat("Extracting water fraction per sounding (exactextractr)...\n")
system.time({
  water_sounding_new <- exactextractr::exact_extract(
    is_water,
    soundings_wgs,
    fun = "mean"
  )
})

water_sounding_new[is.na(water_sounding_new)] <- 0
water_sounding_new <- pmin(pmax(water_sounding_new, 0), 1)

# =============================================================================
# 4. Compare with existing values
# =============================================================================

water_old <- goebel2026::soundings_augmented$proportion_water
water_old[is.na(water_old)] <- 0

cat("\n--- Comparison ---\n")
cat(sprintf("Old:  mean=%.4f  sd=%.4f  prop_zero=%.3f\n",
            mean(water_old), sd(water_old), mean(water_old == 0)))
cat(sprintf("New:  mean=%.4f  sd=%.4f  prop_zero=%.3f\n",
            mean(water_sounding_new), sd(water_sounding_new),
            mean(water_sounding_new == 0)))
cat(sprintf("Correlation (old vs new): %.4f\n",
            cor(water_old, water_sounding_new)))
cat(sprintf("Soundings where old=0 but new>0: %d\n",
            sum(water_old == 0 & water_sounding_new > 0)))

cat("\nQuantiles of new water:\n")
print(quantile(water_sounding_new, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))

cat("\nQuantiles of old water:\n")
print(quantile(water_old, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))

# =============================================================================
# 5. Save
# =============================================================================

out <- data.frame(
  sounding_id   = seq_len(nrow(soundings_proj)),
  water_old     = water_old,
  water_new     = water_sounding_new
)

#write.csv(out, "data-raw/sounding_water_comparison_v7.csv", row.names = FALSE)
#cat("\nSaved to data-raw/sounding_water_comparison_v7.csv\n")

soundings_augmented$proportion_water <- water_sounding_new

# Save just the new water vector for use in modeling scripts
usethis::use_data(soundings_augmented, overwrite = TRUE)
