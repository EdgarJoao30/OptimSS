# This script calculates the focal discrete metrics between the simulated raster and the predicted rasters 
# for mosquito abundance. It classifies the rasters into discrete classes based on quantiles, 
# calculates the frequency of each class within a moving window, and then computes the squared differences 
# between the simulated and predicted frequencies. 
# Focal landscape metrics:
# - class proportion: Proportion of each class within the moving window (composition metric)
# - lsm_l_area_cv: Coefficient of variation of patch area (area and edge metric)
#     - McGarigal K., SA Cushman, and E Ene. 2023. FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical Maps. Computer software program produced by the authors; available at the following web site: https://www.fragstats.org
# - lsm_l_joinent: Join count entropy (complexity metric)
#     - Nowosad J., TF Stepinski. 2019. Information theory as a consistent framework for quantification and classification of landscape patterns. https://doi.org/10.1007/s10980-019-00830-x 
# - lsm_l_contig_mn: Mean patch contiguity (shape metric)
#     - LaGro, J. 1991. Assessing patch shape in landscape mosaics. Photogrammetric Engineering and Remote Sensing, 57(3), 285-293
# - lsm_l_frac_mn: Mean patch fractal dimension (shape metric)
#     - Mandelbrot, B. B. 1977. Fractals: Form, Chance, and Dimension. San Francisco. W. H. Freeman and Company.
# - lsm_l_relmutinf: Relative mutual information (complexity metric)
#     - Nowosad J., TF Stepinski. 2019. Information theory as a consistent framework for quantification and classification of landscape patterns. https://doi.org/10.1007/s10980-019-00830-x
packs <- c("tidyverse", "terra", "sf", "stars", "tidyterra", "patchwork", "landscapemetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

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

calculate_frequency <- function(classified_data, n_classes = 6, w_size = 13, mask_raster = sim) {
  w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
  lapply(0:(n_classes - 1), function(class) {
    dummy_class <- (classified_data == class) * 1
    freq <- focal(dummy_class, w = w_matrix, fun = sum, na.rm = TRUE) / (w_size^2) 
    return(freq |>
             mask(mask_raster[[1]]))
  })
}

calculate_squared_diff <- function(freq_list1, freq_list2) {
  lapply(seq_along(freq_list1), function(idx) {
    (freq_list1[[idx]] - freq_list2[[idx]])^2
  })
}

rmse_focal_discrete_metrics <- function(sim, pred, w_size = 13) {
  # Classify rasters by quantiles
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  
  # Setup parameters
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
  
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
  # avg_diff <- Reduce("+", sum_diff) / n_layers
  # avg_avg_diff <- mean(values(avg_diff), na.rm = TRUE)
  
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
  
  list(
    class_freq = rmse_class_freq,
    joinent = rmse_joinent,
    contig_mn = rmse_contig_mn
  )
}


# Apply function
start_time <- Sys.time()
focal_metrics_a <- rmse_focal_discrete_metrics(sim, a, w_size = 13)
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat("Elapsed time:", elapsed_time, "seconds\n")


focal_metrics_a_27 <- rmse_focal_discrete_metrics(sim, a, w_size = 27)

plot(focal_metrics_a |> rast())
avg_a <- lapply(focal_metrics_a, function(x) mean(values(x), na.rm = TRUE)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "mean_value")

avg_a_27 <- lapply(focal_metrics_a_27, function(x) mean(values(x), na.rm = TRUE)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "mean_value")


focal_metrics_i <- rmse_focal_discrete_metrics(sim, i, w_size = 13)
focal_metrics_i_27 <- rmse_focal_discrete_metrics(sim, i, w_size = 27)

avg_i <- lapply(focal_metrics_i, function(x) mean(values(x), na.rm = TRUE)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "mean_value")

avg_i_27 <- lapply(focal_metrics_i_27, function(x) mean(values(x), na.rm = TRUE)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "mean_value")


plot(focal_metrics_i |> rast())

plot_window_result <- function(data, title, limits, metric_name = "relmutinf") {
  ggplot() +
    geom_spatraster(data = data) +
    scale_fill_viridis_c(limits = limits, oob = scales::squish, na.value = NA) +
    labs(title = title, fill = paste0("RMSE ", metric_name)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    )
}

# Determine the common color scale limits
common_limits <- range(c(values(focal_metrics_a$relmutinf), values(focal_metrics_i$relmutinf)), na.rm = TRUE)

# Create the plots
plot_a <- plot_window_result(focal_metrics_a$relmutinf, "Uniform-Fixed\nWindow size: 13", common_limits, metric_name = "relmutinf")
plot_i <- plot_window_result(focal_metrics_i$relmutinf, "Random-Variable\nWindow size: 13", common_limits, metric_name = "relmutinf")

# Combine the plots side by side
combined_plot <- plot_a + plot_i + plot_layout(ncol = 2, guides = 'collect')

# Display the combined plot
print(combined_plot)
