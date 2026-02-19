packs <- c("tidyverse", "terra", "sf", "stars")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
# Function to classify continuous rasters into categorical rasters based on quantiles
classify_raster_by_quantiles <- function(raster_data, n_breaks = 6) {
  if (!inherits(raster_data, "SpatRaster")) {
    stop("raster_data must be a SpatRaster.")
  }

  n_layers <- nlyr(raster_data)
  classified_layers <- lapply(seq_len(n_layers), function(layer_idx) {
    layer <- raster_data[[layer_idx]]
    values_vec <- layer[!is.na(layer)]
    breaks <- quantile(values_vec, probs = seq(0, 1, length.out = n_breaks), na.rm = TRUE)
    class_matrix <- matrix(
      c(breaks[-length(breaks)], breaks[-1], seq_along(breaks[-1]) - 1),
      ncol = 3,
      byrow = FALSE
    )
    class_matrix[1, 1] <- 0
    classify(layer, class_matrix)
  })

  rast(classified_layers)
}
# Function to calculate focal class frequency for a given window size
calculate_frequency <- function(classified_data, n_classes = 6, w_size = 13, mask_raster = sim) {
  w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
  lapply(0:(n_classes - 1), function(class) {
    dummy_class <- (classified_data == class) * 1
    freq <- focal(dummy_class, w = w_matrix, fun = sum, na.rm = TRUE) / (w_size^2) 
    return(freq |>
             mask(mask_raster[[1]]))
  })
}
# Function to calculate squared differences between two raster time series
calculate_squared_diff <- function(freq_list1, freq_list2) {
  lapply(seq_along(freq_list1), function(idx) {
    (freq_list1[[idx]] - freq_list2[[idx]])^2
  })
}
# Function to calculate kolmogorov-smirnov distances on patch level metrics
# inputs are the outputs of calculate_lsm for the simulated and predicted rasters, and the name of the metric to compare
calculate_ks_metric <- function(sim_metrics, pred_metrics, metric_name, n_layers) {
  sim_metric <- sim_metrics |> filter(metric == metric_name)
  pred_metric <- pred_metrics |> filter(metric == metric_name)
  ks_results <- lapply(1:n_layers, function(i) {
    sim_values <- sim_metric |> filter(layer == i) |> pull(value)
    pred_values <- pred_metric |> filter(layer == i) |> pull(value)
    ks_result <- ks.test(sim_values, pred_values)
    return(ks_result$statistic)
  })
  return(ks_results |> unlist() |> mean())
}