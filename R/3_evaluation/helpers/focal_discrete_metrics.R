# This script calculates the focal discrete metrics between the simulated raster and the predicted rasters 
# for mosquito abundance. It classifies the rasters into discrete classes based on quantiles, 
# calculates the frequency of each class within a moving window, and then computes the squared differences 
# between the simulated and predicted frequencies. 
packs <- c("tidyverse", "terra", "sf", "stars", "tidyterra", "patchwork", "landscapemetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

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

# Apply classification function
sim_classified <- classify_raster_by_quantiles(sim)
a_classified <- classify_raster_by_quantiles(a)
i_classified <- classify_raster_by_quantiles(i)


n_classes <- length((unique(values(sim_classified[[1]], na.rm = TRUE))))
n_layers <- nlyr(sim_classified)
w_size <- 15
w_matrix <- matrix(1, nrow = w_size, ncol = w_size)

calculate_frequency <- function(classified_data) {
  lapply(0:(n_classes - 1), function(class) {
    dummy_class <- (classified_data == class) * 1
    freq <- focal(dummy_class, w = w_matrix, fun = sum, na.rm = TRUE) / (w_size^2) 
    return(freq |>
             mask(sim[[1]]))
  })
}

sim_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(sim_classified[[i]]))
a_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(a_classified[[i]]))
i_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(i_classified[[i]]))

calculate_squared_diff <- function(freq_list1, freq_list2) {
  lapply(seq_along(freq_list1), function(idx) {
    (freq_list1[[idx]] - freq_list2[[idx]])^2
  })
}

diff_a <- lapply(seq_len(n_layers), function(i) calculate_squared_diff(sim_freq_list[[i]], a_freq_list[[i]]))
diff_i <- lapply(seq_len(n_layers), function(i) calculate_squared_diff(sim_freq_list[[i]], i_freq_list[[i]]))

calculate_per_class_diff <- function(diff_list) {
  lapply(seq_len(n_layers), function(i) {
    sapply(diff_list[[i]], function(layer) {
      mean(values(layer), na.rm = TRUE)
    })
  })
}

per_class_a <- calculate_per_class_diff(diff_a)
per_class_i <- calculate_per_class_diff(diff_i)

per_class_a_mean <- sapply(seq_len(n_classes), function(j) {
  mean(sapply(seq_len(n_layers), function(i) per_class_a[[i]][j]))
}) 

per_class_i_mean <- sapply(seq_len(n_classes), function(j) {
  mean(sapply(seq_len(n_layers), function(i) per_class_i[[i]][j]))
}) 

sum_a <- lapply(seq_len(n_layers), function(i) Reduce("+", diff_a[[i]]))
sum_i <- lapply(seq_len(n_layers), function(i) Reduce("+", diff_i[[i]]))

avg_a <- Reduce("+", sum_a) / n_layers
avg_i <- Reduce("+", sum_i) / n_layers

avg_avg_a <- mean(values(avg_a), na.rm = TRUE)
avg_avg_i <- mean(values(avg_i), na.rm = TRUE)


calculate_focal_discrete_metrics <- function(sim, pred, w_size = 15) {
  # Classify rasters by quantiles
  sim_classified <- classify_raster_by_quantiles(sim)
  pred_classified <- classify_raster_by_quantiles(pred)
  
  # Setup parameters
  n_classes <- length(unique(values(sim_classified[[1]], na.rm = TRUE)))
  n_layers <- nlyr(sim_classified)
  w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
  
  # Calculate frequency for each class within moving window
  sim_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(sim_classified[[i]]))
  pred_freq_list <- lapply(seq_len(n_layers), function(i) calculate_frequency(pred_classified[[i]]))
  
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
  
  # Calculate cumulative and average differences
  sum_diff <- lapply(seq_len(n_layers), function(i) Reduce("+", diff_list[[i]]))
  avg_diff <- Reduce("+", sum_diff) / n_layers
  avg_avg_diff <- mean(values(avg_diff), na.rm = TRUE)
  
  list(
    diff = diff_list,
    per_class_mean = per_class_mean,
    sum_diff = sum_diff,
    avg_diff = avg_diff,
    avg_avg_diff = avg_avg_diff
  )
}

# Apply function
focal_metrics_a <- calculate_focal_discrete_metrics(sim, a, w_size = 15)
focal_metrics_i <- calculate_focal_discrete_metrics(sim, i, w_size = 15)
diff_a <- focal_metrics_a$diff
plot(diff_a[[1]][[1]])
per_class_a_mean <- focal_metrics_a$per_class_mean
sum_a <- focal_metrics_a$sum_diff
avg_a <- focal_metrics_a$avg_diff
plot(avg_a)
avg_avg_a <- focal_metrics_a$avg_avg_diff

diff_i <- focal_metrics_i$diff
per_class_i_mean <- focal_metrics_i$per_class_mean
sum_i <- focal_metrics_i$sum_diff
avg_i <- focal_metrics_i$avg_diff
plot(avg_i)
avg_avg_i <- focal_metrics_i$avg_avg_diff

plot_window_result <- function(data, title, limits) {
  ggplot() +
    geom_spatraster(data = data) +
    scale_fill_viridis_c(limits = limits, oob = scales::squish, na.value = NA) +
    labs(title = title, fill = "Distance") +
    theme_void(base_family = 'Times New Roman', base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    )
}

# Determine the common color scale limits
common_limits <- range(c(values(sum_a), values(sum_i)), na.rm = TRUE)

# Create the plots
plot_a <- plot_window_result(sum_a, "Window Result A", common_limits)
plot_i <- plot_window_result(sum_i, "Window Result I", common_limits)

# Combine the plots side by side
combined_plot <- plot_a + plot_i + plot_layout(ncol = 2, guides = 'collect')

# Display the combined plot
print(combined_plot)
