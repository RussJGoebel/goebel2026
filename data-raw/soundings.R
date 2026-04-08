## code to prepare `soundings_synthetic_albedo` dataset goes here

### Setup ######################################################################
set.seed(1) # the upscaling requires random noise
library(spatintegrate)
library(downscaling)
soundings_synthetic_albedo <- downscaling::soundings # the soundings in the paper are the same as in downscaling package
soundings_synthetic_albedo <- downscaling::soundings #ensure_projected(downscaling::soundings)

future::plan(future::multisession())


### /Setup #####################################################################

### Parameters #################################################################


### /Parameters ################################################################

### Upscaling ##################################################################

na_indices <- (is.na(target_grid$mean_albedo))
noise_sd <- sd(target_grid$mean_albedo[!na_indices])/20  #we set the standard deviation of the noise to about 5% of the population

A <- spatintegrate::compute_overlap_fractions(ensure_projected(soundings),
                                              ensure_projected(target_grid))

mean_value <- target_grid$mean_albedo[!na_indices]
e <- rnorm(dim(A)[1],0,noise_sd)

soundings_synthetic_albedo$synthetic_albedo_no_noise <- as.vector(A[,!na_indices] %*% mean_value)

### /Upscaling #################################################################

# -----------------------------------------------

### Save data:

usethis::use_data(soundings_synthetic_albedo, overwrite = TRUE)


###

# median_albedo_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_albedo_june_2022_largebbox.tif")
# median_ndvi_rast <- terra::rast("/projectnb/buultra/SIF_downscaling/russell/data/median_boston_ndvi_june_2022.tif")
#
# soundings_synthetic_albedo2 <- soundings_synthetic_albedo
#
# soundings_synthetic_albedo2  <- summarize_raster_mean_over_polygons(median_ndvi_rast,soundings_synthetic_albedo2)
#
# soundings_synthetic_albedo2  <- summarize_raster_mean_over_polygons(median_albedo_rast,soundings_synthetic_albedo2)
# soundings_synthetic_albedo2$mean_albedo <- soundings_synthetic_albedo2$mean_value
# soundings_synthetic_albedo2$mean_value <- NULL
#
# soundings_synthetic_albedo <- soundings_synthetic_albedo2
#
# usethis::use_data(soundings_synthetic_albedo, overwrite = TRUE)
