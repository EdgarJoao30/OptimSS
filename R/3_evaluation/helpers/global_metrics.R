packs <- c("tidyverse", "terra", "sf", "stars", "landscapemetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")
a <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/a_15_allmonths.geojson") |>
  dplyr::filter(variable == "pred_mean") |>
  dplyr::select(month, value, geometry) |>
  split(~month) |>
  lapply(function(x) {
    x |>
      dplyr::select(value, geometry) |>
      st_rasterize() |>
      rast() |>
      mask(sim[[1]])
  }) |>
  rast()  

pred <- a
global_metrics <- function(sim, pred) {
  ### Settings
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  metrics <- c("lsm_l_joinent", "lsm_l_contig_mn",
               "lsm_p_area", "lsm_p_enn", "lsm_p_frac")
  # Abundance RMSE

  # Class frequency differences
  
  # Landscape and patch metrics
  sim_metrics <- calculate_lsm(sim_classified, what = metrics)
  pred_metrics <- calculate_lsm(pred_classified, what = metrics)

  sim_metrics
  
  return(list())
}