library(terra)
library(sf)
library(stars)
library(tidyterra)
library(tidyverse)
library(patchwork)

sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")

a <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/a_15_an_df.geojson") |>
  dplyr::filter(month == 1 & variable == "pred_mean") |>
  dplyr::select(value, geometry) |>
  st_rasterize() |>
  rast() |>
  mask(sim[[1]])

i <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/i_50_an_df.geojson") |>
  dplyr::filter(month == 1 & variable == "pred_mean") |>
  dplyr::select(value, geometry) |>
  st_rasterize() |>
  rast() |>
  mask(sim[[1]])

focal_euclidean <- function(r1, r2, window_size = 7) {

  # Check if rasters align
  compareGeom(r1, r2)

  # Step A: Calculate Squared Difference (Pixel-by-Pixel)
  # (Map1 - Map2)^2
  diff_sq <- (r1 - r2)^2

  # Step B: Define the Window (Weights Matrix)
  # A matrix of 1s indicates we sum all pixels in the window equally
  w_matrix <- matrix(1, nrow=window_size, ncol=window_size)

  # Step C: Apply Focal Sum
  # This sums the squared differences within the window neighborhood
  # na.rm=TRUE handles edges or missing data gracefully
  sum_sq_window <- focal(diff_sq, w=w_matrix, fun=sum, na.rm=TRUE)

  # Step D: Take the Square Root to get Euclidean Distance
  local_dist <- sqrt(sum_sq_window)
  return(local_dist)
}

window_result_a <- focal_euclidean(sim[[1]], a, window_size = 14)
window_result_i <- focal_euclidean(sim[[1]], i, window_size = 14)
plot(window_result_a)
plot(window_result_i)

# Define a common plotting function for the window results
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
common_limits <- range(c(values(window_result_a), values(window_result_i)), na.rm = TRUE)

# Create the plots
plot_a <- plot_window_result(window_result_a, "Window Result A", common_limits)
plot_i <- plot_window_result(window_result_i, "Window Result I", common_limits)

# Combine the plots side by side
combined_plot <- plot_a + plot_i + plot_layout(ncol = 2, guides = 'collect')

# Display the combined plot
print(combined_plot)
