## =============================================================================
## prepare_target_grid.R
##
## Builds target_grid: a 330m regular lattice over Boston augmented with:
##   - mean_albedo       (Sentinel-2 median albedo, June 2022)
##   - mean_ndvi         (Sentinel-2 median NDVI, June 2022)
##   - landcover class proportions (ESA WorldCover via existing spatintegrate fns)
##   - n_intersects      (number of soundings overlapping each cell)
##
## Uses exactextractr for albedo and NDVI extraction.
## Requires: module load geos before starting R
## =============================================================================

library(goebel2026)
library(spatintegrate)
library(terra)
library(sf)
library(exactextractr)

### Parameters #################################################################

resolution_of_grid <- 330
buffer             <- 25 * resolution_of_grid

landcover_key <- c(
  "1" = "urban",
  "2" = "croplands",
  "3" = "forest",
  "4" = "shrub & scrub",
  "5" = "grass",
  "6" = "water",
  "7" = "wetlands",
  "8" = "bare surface"
)

### Load rasters ###############################################################

landcover_rast  <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif")
albedo_rast     <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif")
ndvi_rast       <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_ndvi_june_2022.tif")

### 1) Create grid #############################################################

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)

target_grid <- spatintegrate::make_square_grid_in_crs(
  soundings_proj,
  resolution_of_grid,
  buffer = buffer
)

# Transform to raster CRS for extraction, then we'll reproject back at the end
target_grid_wgs <- sf::st_transform(target_grid, sf::st_crs(goebel2026::soundings))

### 2) Crop rasters to grid extent ############################################

cat("Cropping rasters to target grid extent...\n")
grid_vect_lc     <- terra::vect(sf::st_transform(target_grid_wgs, sf::st_crs(landcover_rast)))
grid_vect_albedo <- terra::vect(sf::st_transform(target_grid_wgs, sf::st_crs(albedo_rast)))
grid_vect_ndvi   <- terra::vect(sf::st_transform(target_grid_wgs, sf::st_crs(ndvi_rast)))

landcover_crop <- terra::crop(landcover_rast, terra::ext(grid_vect_lc))
albedo_crop    <- terra::crop(albedo_rast,    terra::ext(grid_vect_albedo))
ndvi_crop      <- terra::crop(ndvi_rast,      terra::ext(grid_vect_ndvi))

cat(sprintf("Landcover crop:  %d x %d\n", nrow(landcover_crop), ncol(landcover_crop)))
cat(sprintf("Albedo crop:     %d x %d\n", nrow(albedo_crop),    ncol(albedo_crop)))
cat(sprintf("NDVI crop:       %d x %d\n", nrow(ndvi_crop),      ncol(ndvi_crop)))

### 3) Extract albedo and NDVI via exactextractr ###############################

cat("Extracting mean_albedo...\n")
grid_albedo <- sf::st_transform(target_grid_wgs, sf::st_crs(albedo_rast))
system.time({
  target_grid_wgs$mean_albedo <- exactextractr::exact_extract(albedo_crop, grid_albedo, fun = "mean")
})

cat("Extracting mean_ndvi...\n")
grid_ndvi <- sf::st_transform(target_grid_wgs, sf::st_crs(ndvi_rast))
system.time({
  target_grid_wgs$mean_ndvi <- exactextractr::exact_extract(ndvi_crop, grid_ndvi, fun = "mean")
})

cat(sprintf("mean_albedo: mean=%.4f  NAs=%d\n", mean(target_grid_wgs$mean_albedo, na.rm=TRUE), sum(is.na(target_grid_wgs$mean_albedo))))
cat(sprintf("mean_ndvi:   mean=%.4f  NAs=%d\n", mean(target_grid_wgs$mean_ndvi,   na.rm=TRUE), sum(is.na(target_grid_wgs$mean_ndvi))))

### 4) Augment grid with landcover ############################################

# These functions do class-proportion extraction -- keeping as-is
grid_lc <- sf::st_transform(target_grid_wgs, sf::st_crs(landcover_rast))

target_grid_wgs <- spatintegrate::summarize_raster_class_representation_over_grid(
  landcover_crop,
  target_grid_wgs
)
target_grid_wgs <- goebel2026::rename_and_fill_proportions(
  target_grid_wgs,
  key = landcover_key
)
target_grid_wgs <- goebel2026::rename_and_fill_proportions(
  target_grid_wgs,
  prefix = "n_",
  key = landcover_key
)

### 5) Count sounding intersections per grid cell #############################

target_grid_wgs$n_intersects <- lengths(
  sf::st_intersects(target_grid_wgs, goebel2026::soundings)
)

### 6) Reproject back to soundings CRS and save ################################

target_grid <- sf::st_transform(
  target_grid_wgs,
  crs = sf::st_crs(goebel2026::soundings)
)

usethis::use_data(target_grid, overwrite = TRUE)
cat("Saved target_grid.\n")
