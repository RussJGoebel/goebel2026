# make_target_grid_square.R
#
# Generates target_grid_square: a square grid centered on the soundings centroid.
# Side length L = max(Lx, Ly) of the original rectangular grid, extended to
# a square so that SAR eigenvectors exactly match the Fourier basis on [0,L]^2.
#
# The extra pixels (outside the original rectangle) are unobserved and
# prior-dominated -- they don't affect the fit inside the original rectangle.

library(goebel2026)
library(spatintegrate)
library(sf)

### Parameters #################################################################
resolution_of_grid <- 330
buffer <- 25 * resolution_of_grid

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

### 1) Create rectangular grid (same as original) #############################
soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)

target_grid_rect <- spatintegrate::make_square_grid_in_crs(
  soundings_proj,
  resolution_of_grid,
  buffer = buffer
)

# Get rectangular bbox and dimensions
rect_bbox <- as.numeric(sf::st_bbox(target_grid_rect))
Lx <- rect_bbox[3] - rect_bbox[1]
Ly <- rect_bbox[4] - rect_bbox[2]
m1 <- round(Ly / resolution_of_grid)
m2 <- round(Lx / resolution_of_grid)

cat(sprintf("Original grid: m1=%d rows x m2=%d cols  (Lx=%.0f, Ly=%.0f)\n",
            m1, m2, Lx, Ly))

### 2) Build square grid CENTERED on soundings centroid #######################
J <- max(m1, m2)
L <- J * resolution_of_grid

cat(sprintf("Square grid:   J=%d x %d  (L=%.0f)\n", J, J, L))

# Centroid of soundings in projected CRS
ctr <- sf::st_centroid(sf::st_union(sf::st_geometry(soundings_proj))) |>
  sf::st_coordinates()

cat(sprintf("Soundings centroid: x=%.0f  y=%.0f\n", ctr[1], ctr[2]))

# Square bbox centered on soundings centroid
sq_bbox <- c(
  xmin = ctr[1] - L/2,
  ymin = ctr[2] - L/2,
  xmax = ctr[1] + L/2,
  ymax = ctr[2] + L/2
)

cat(sprintf("Square bbox: xmin=%.0f ymin=%.0f xmax=%.0f ymax=%.0f\n",
            sq_bbox["xmin"], sq_bbox["ymin"],
            sq_bbox["xmax"], sq_bbox["ymax"]))

# Check soundings are inside square
snd_bbox <- as.numeric(sf::st_bbox(soundings_proj))
cat(sprintf("Soundings bbox: xmin=%.0f ymin=%.0f xmax=%.0f ymax=%.0f\n",
            snd_bbox[1], snd_bbox[2], snd_bbox[3], snd_bbox[4]))
cat(sprintf("Soundings inside square: %s\n",
            all(snd_bbox[1] >= sq_bbox["xmin"],
                snd_bbox[2] >= sq_bbox["ymin"],
                snd_bbox[3] <= sq_bbox["xmax"],
                snd_bbox[4] <= sq_bbox["ymax"])))

# Generate square grid
sq_grid_sfc <- sf::st_make_grid(
  sf::st_as_sfc(sf::st_bbox(sq_bbox, crs = sf::st_crs(soundings_proj))),
  cellsize = resolution_of_grid,
  what     = "polygons",
  square   = TRUE
)
target_grid_square <- sf::st_sf(geometry = sq_grid_sfc)

cat(sprintf("Square grid cells: %d  (expected J^2=%d)\n",
            nrow(target_grid_square), J^2))

### 3) Transform for raster extraction ########################################
target_grid_sq_transform <- sf::st_transform(
  target_grid_square,
  sf::st_crs(goebel2026::soundings)
)

### 4) Augment with median albedo #############################################
target_grid_sq_transform <- spatintegrate::summarize_raster_mean_over_polygons(
  median_albedo_rast,
  target_grid_sq_transform,
  stats_col_prefix = "mean_"
)
names(target_grid_sq_transform)[
  names(target_grid_sq_transform) == "mean_albedo_median"] <- "mean_albedo"

### 5) Augment with median NDVI ###############################################
target_grid_sq_transform <- spatintegrate::summarize_raster_mean_over_polygons(
  median_ndvi_rast,
  target_grid_sq_transform,
  stats_col_prefix = "mean_"
)
names(target_grid_sq_transform)[
  names(target_grid_sq_transform) == "mean_NDVI_median"] <- "mean_ndvi"

### 6) Augment with landcover #################################################
target_grid_sq_transform <- spatintegrate::summarize_raster_class_representation_over_grid(
  landcover_rast,
  target_grid_sq_transform
)
target_grid_sq_transform <- goebel2026::rename_and_fill_proportions(
  target_grid_sq_transform,
  key = landcover_key
)
target_grid_sq_transform <- goebel2026::rename_and_fill_proportions(
  target_grid_sq_transform,
  prefix = "n_",
  key = landcover_key
)

### 7) Count sounding intersections ###########################################
target_grid_sq_transform$n_intersects <- lengths(
  sf::st_intersects(target_grid_sq_transform, goebel2026::soundings)
)

# Flag which cells are in original rectangle vs padding
rect_bbox_sfc <- sf::st_transform(
  sf::st_as_sfc(sf::st_bbox(target_grid_rect)),
  sf::st_crs(goebel2026::soundings)
)
target_grid_sq_transform$in_original <- as.numeric(
  sf::st_intersects(target_grid_sq_transform, rect_bbox_sfc)
) > 0

cat(sprintf("Cells in original rectangle: %d\n",
            sum(target_grid_sq_transform$in_original, na.rm=TRUE)))
cat(sprintf("Padding cells:               %d\n",
            sum(!target_grid_sq_transform$in_original, na.rm=TRUE)))

### 8) Reproject back #########################################################
target_grid_square <- sf::st_transform(
  target_grid_sq_transform,
  crs = sf::st_crs(goebel2026::soundings)
)

### 9) Verify square and centring #############################################
sq_bbox_final <- as.numeric(sf::st_bbox(
  spatintegrate::ensure_projected(target_grid_square)
))
Lx_sq <- sq_bbox_final[3] - sq_bbox_final[1]
Ly_sq <- sq_bbox_final[4] - sq_bbox_final[2]
ctr_sq <- c((sq_bbox_final[1]+sq_bbox_final[3])/2,
            (sq_bbox_final[2]+sq_bbox_final[4])/2)

cat(sprintf("\nFinal square grid: Lx=%.0f Ly=%.0f (should be equal)\n",
            Lx_sq, Ly_sq))
cat(sprintf("m1_sq=%.0f m2_sq=%.0f (should both = %d)\n",
            Ly_sq/resolution_of_grid, Lx_sq/resolution_of_grid, J))
cat(sprintf("Grid centroid:     x=%.0f  y=%.0f\n", ctr_sq[1], ctr_sq[2]))
cat(sprintf("Soundings centroid: x=%.0f  y=%.0f\n", ctr[1], ctr[2]))
cat(sprintf("Offset: dx=%.0f  dy=%.0f  (should be ~0)\n",
            ctr_sq[1]-ctr[1], ctr_sq[2]-ctr[2]))

### Save ######################################################################
usethis::use_data(target_grid_square, overwrite = TRUE)
cat("Saved target_grid_square\n")
