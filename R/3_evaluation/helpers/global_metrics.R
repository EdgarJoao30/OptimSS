packs <- c("tidyverse", "terra", "sf", "stars", "landscapemetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

# sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")
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

# pred <- a
global_metrics <- function(sim, pred) {
  ### Settings
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  metrics <- c("lsm_l_joinent", "lsm_l_contig_mn",
               "lsm_p_area", "lsm_p_enn", "lsm_p_frac")
  # Abundance RMSE
  sim_mean <- sim |> terra::global(mean, na.rm = TRUE) |> pull(mean)
  pred_mean <- pred |> terra::global(mean, na.rm = TRUE) |> pull(mean)
  abundance_rmse <- sqrt(mean((sim_mean - pred_mean)^2))
  
  # Landscape and patch metrics
  sim_metrics <- calculate_lsm(sim_classified, what = metrics)
  pred_metrics <- calculate_lsm(pred_classified, what = metrics)

  sim_contig_mn <- sim_metrics |> filter(metric == "contig_mn") |> pull(value)
  pred_contig_mn <- pred_metrics |> filter(metric == "contig_mn") |> pull(value)
  contig_rmse <- sqrt(mean((pred_contig_mn - sim_contig_mn)^2))
  
  sim_joinent <- sim_metrics |> filter(metric == "joinent") |> pull(value)
  pred_joinent <- pred_metrics |> filter(metric == "joinent") |> pull(value)
  joinent_rmse <- sqrt(mean((pred_joinent - sim_joinent)^2))
  
  ks_area <- calculate_ks_metric(sim_metrics = sim_metrics, pred_metrics = pred_metrics, metric_name = "area", n_layers = n_layers)
  ks_enn <- calculate_ks_metric(sim_metrics = sim_metrics, pred_metrics = pred_metrics, metric_name = "enn", n_layers = n_layers)
  ks_frac <- calculate_ks_metric(sim_metrics = sim_metrics, pred_metrics = pred_metrics, metric_name = "frac", n_layers = n_layers)

  return(list(global_abundance_rmse = abundance_rmse,
              global_joinent_rmse = joinent_rmse,
              global_contig_rmse = contig_rmse,
              global_area_ks = ks_area,
              global_enn_ks = ks_enn,
              global_frac_ks = ks_frac))
}

# test
# test <- global_metrics(sim, pred)
