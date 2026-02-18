# This script calculates the focal discrete metrics between the simulated raster and the predicted rasters 
# for mosquito abundance. It classifies the rasters into discrete classes based on quantiles, 
# calculates the frequency of each class within a moving window, and then computes the squared differences 
# between the simulated and predicted frequencies. 
# Focal abundance RMSE is also calculated by applying a moving window to the original rasters and 
# computing the RMSE between the windowed means.
# Focal landscape metrics:
# - class proportion: Proportion of each class within the moving window (composition metric)
# - lsm_l_joinent: Join count entropy (complexity metric)
#     - Nowosad J., TF Stepinski. 2019. Information theory as a consistent framework for quantification and classification of landscape patterns. https://doi.org/10.1007/s10980-019-00830-x 
# - lsm_l_contig_mn: Mean patch contiguity (shape metric)
#     - LaGro, J. 1991. Assessing patch shape in landscape mosaics. Photogrammetric Engineering and Remote Sensing, 57(3), 285-293
packs <- c("tidyverse", "terra", "sf", "stars", "tidyterra", "patchwork", "landscapemetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")
mean_sim <- rast("~/Documents/GitHub/OptimSS/data/1_raw/20260214_sim_raster001_mean.tif")
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

i <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/i_15_allmonths.geojson") |>
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
focal_metrics <- function(sim, pred, w_size = 13) {
  ### Settings
  # Classify rasters by quantiles
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  w_matrix <- matrix(1, nrow = w_size, ncol = w_size)

  ### Continuous metrics
  # Calculate squared differences
  mean_sim <- focal(sim, w = w_matrix, fun = mean, na.rm = TRUE) 
  mean_pred <- focal(pred, w = w_matrix, fun = mean, na.rm = TRUE) 
  diff_sq <- (mean_sim - mean_pred)^2
  rmse_window <- sqrt(mean(diff_sq, na.rm = TRUE))
  rmse_window <- rmse_window |> mask(sim[[1]])
  overall_rmse <- mean(values(rmse_window), na.rm = TRUE)
  ### Discrete metrics
  # Calculate frequency for each class within moving window
  sim_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(sim_classified[[i]], n_classes, w_size, mask_raster = sim[[1]]))
  pred_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(pred_classified[[i]], n_classes, w_size, mask_raster = sim[[1]]))
  
  # Calculate squared differences
  diff_list <- lapply(seq_len(n_layers), function(i) calculate_squared_diff(sim_freq_list[[i]], pred_freq_list[[i]]))
  
  # Calculate per-class differences
  per_class_diff <- lapply(seq_len(n_layers), function(i) {
    sapply(diff_list[[i]], function(layer) {
      mean(values(layer), na.rm = TRUE)
    })
  })
  
  per_class_mean <- sapply(seq_len(n_classes), function(j) {
    mean(sapply(seq_len(n_layers), function(i) per_class_diff[[i]][j]))
  })
  
  # Calculate class frequency cumulative differences
  class_freq <- lapply(seq_len(n_layers), function(i) Reduce("+", diff_list[[i]]))
  
  # Calculate landscape metrics 
  sim_lm <- window_lsm(sim_classified, window = w_matrix, 
                       what = c("lsm_l_joinent", "lsm_l_contig_mn"))
  pred_lm <- window_lsm(pred_classified, window = w_matrix, 
                        what = c("lsm_l_joinent", "lsm_l_contig_mn"))
  
  metrics <- c("lsm_l_joinent", "lsm_l_contig_mn")

  for (metric in metrics) {
    assign(paste0("sim_", gsub("lsm_l_", "", metric)), 
           lapply(sim_lm, function(x) x[[metric]] |> mask(sim[[1]])))
  }
  
  for (metric in metrics) {
    assign(paste0("pred_", gsub("lsm_l_", "", metric)), 
           lapply(pred_lm, function(x) x[[metric]] |> mask(sim[[1]])))
  }
  
  diff_joinent <- calculate_squared_diff(sim_joinent, pred_joinent) 
  diff_contig_mn <- calculate_squared_diff(sim_contig_mn, pred_contig_mn) 
  
  rmse_class_freq <- sqrt(Reduce("+", class_freq) / n_layers)
  rmse_joinent <- sqrt(Reduce("+", diff_joinent) / n_layers)
  rmse_contig_mn <- sqrt(Reduce("+", diff_contig_mn) / n_layers)
  
  overall_rmse_class_freq <- mean(values(rmse_class_freq), na.rm = TRUE)
  overall_rmse_joinent <- mean(values(rmse_joinent), na.rm = TRUE)
  overall_rmse_contig_mn <- mean(values(rmse_contig_mn), na.rm = TRUE)
  list(
    overall_rmse_window = overall_rmse,
    overall_class_freq = overall_rmse_class_freq,
    overall_joinent = overall_rmse_joinent,
    overall_contig_mn = overall_rmse_contig_mn
  )
}


# Apply function
start_time <- Sys.time()
focal_metrics_i <- focal_metrics(sim, i, w_size = 27)
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat("Elapsed time:", elapsed_time, "seconds\n")