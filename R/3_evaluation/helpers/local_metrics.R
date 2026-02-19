# This script defines a function to calculate local metrics for evaluating the performance 
# of spatial predictions of mosquito abundance. The function takes in two SpatRaster objects 
# (simulated and predicted) and computes the following metrics:
# - Local RMSE: The root mean squared error calculated at each pixel location, 
# averaged across the entire raster.
# - Kappa coefficient: A measure of agreement between the classified simulated and predicted 
# rasters, calculated for each layer and averaged across layers.
packs <- c("tidyverse", "terra", "sf", "stars", "irr")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

# sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")
# mean_sim <- rast("~/Documents/GitHub/OptimSS/data/1_raw/20260214_sim_raster001_mean.tif")
# a <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/a_15_allmonths.geojson") |>
#   dplyr::filter(variable == "pred_mean") |>
#   dplyr::select(month, value, geometry) |>
#   split(~month) |>
#   lapply(function(x) {
#     x |>
#       dplyr::select(value, geometry) |>
#       st_rasterize() |>
#       rast() |>
#       mask(sim[[1]])
#   }) |>
#   rast()  

# i <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/i_15_allmonths.geojson") |>
#   dplyr::filter(variable == "pred_mean") |>
#   dplyr::select(month, value, geometry) |>
#   split(~month) |>
#   lapply(function(x) {
#     x |>
#       dplyr::select(value, geometry) |>
#       st_rasterize() |>
#       rast() |>
#       mask(sim[[1]])
#   }) |>
#   rast()

local_metrics <- function(sim, pred) { 
  ### Settings
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  # Abundance RMSE
  diff_sq <- (sim - pred)^2
  local_rmse <- sqrt(mean(diff_sq, na.rm = TRUE))
  overall_rmse <- mean(values(local_rmse), na.rm = TRUE)

  # Calculate kappa coefficient
  kappa_coeffs <- lapply(seq_len(n_layers), function(i) {
    sim_values <- values(sim_classified[[i]], na.rm = TRUE)
    pred_values <- values(pred_classified[[i]], na.rm = TRUE)
    kappa <- irr::kappa2(data.frame(sim_values, pred_values))$value
    return(kappa)
  })
  overall_kappa <- mean(unlist(kappa_coeffs), na.rm = TRUE)
  return(list(local_abundance_rmse = overall_rmse, local_kappa = overall_kappa))
}
