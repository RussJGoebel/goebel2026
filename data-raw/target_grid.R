## code to prepare `target_grid` dataset goes here

### target_grid consists of a regular lattice located over Boston.
# this code prepares target_grid by:
# 1) generating a grid of cells using Alber's Equal Area projection centered at
# the centroid of the union of the soundings of interest in Boston. (the
# 'soundings' dataset.)
# 2) augmenting the grid using median albedo values gathered from Google Earth Engine
# 3) augmenting the grid using median NDVI values gathered from Google Earth Engine
# 4) augmenting the grid using landcover values gathered from Google Earth Engine,
# specifically the I95 WorldCover data.
###############################################################################

library(goebel2026)
library(spatintegrate)

### Parameters #################################################################

resolution_of_grid <- 330 # grid cell size (in meters)
buffer <- 25 * resolution_of_grid # buffer to add to edges of grid

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

### Load raster data ###########################################################

landcover_rast     <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/WorldCover_I95_large_bbox.tif")
median_albedo_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif")
median_ndvi_rast   <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_ndvi_june_2022.tif")

### 1) Create grid #############################################################

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)

target_grid <- spatintegrate::make_square_grid_in_crs(
  soundings_proj,
  resolution_of_grid,
  buffer = buffer
)

# Transform to soundings CRS for raster extraction
target_grid_transform <- sf::st_transform(target_grid, sf::st_crs(goebel2026::soundings))

### 2) Augment grid with median albedo #########################################

target_grid_transform <- spatintegrate::summarize_raster_mean_over_polygons(
  median_albedo_rast,
  target_grid_transform,
  stats_col_prefix = "mean_"
)
names(target_grid_transform)[names(target_grid_transform) == "mean_albedo_median"] <- "mean_albedo"

### 3) Augment grid with median NDVI ##########################################

target_grid_transform <- spatintegrate::summarize_raster_mean_over_polygons(
  median_ndvi_rast,
  target_grid_transform,
  stats_col_prefix = "mean_"
)
names(target_grid_transform)[names(target_grid_transform) == "mean_NDVI_median"] <- "mean_ndvi"

### 4) Augment grid with landcover ############################################

target_grid_transform <- spatintegrate::summarize_raster_class_representation_over_grid(
  landcover_rast,
  target_grid_transform
)

target_grid_transform <- goebel2026::rename_and_fill_proportions(
  target_grid_transform,
  key = landcover_key
)
target_grid_transform <- goebel2026::rename_and_fill_proportions(
  target_grid_transform,
  prefix = "n_",
  key = landcover_key
)

### 5) Count sounding intersections per grid cell #############################

target_grid_transform$n_intersects <- lengths(
  sf::st_intersects(target_grid_transform, goebel2026::soundings)
)

### 6) Reproject back to soundings CRS ########################################

target_grid <- sf::st_transform(
  target_grid_transform,
  crs = sf::st_crs(goebel2026::soundings)
)

### Save ######################################################################
usethis::use_data(target_grid, overwrite = TRUE)
