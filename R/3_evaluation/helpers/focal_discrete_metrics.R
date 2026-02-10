library(terra)
library(sf)
library(stars)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(landscapemetrics)

sim <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/20250331_sim_raster001.tif")

a <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/a_15_an_df.geojson") |>
  dplyr::filter(variable == "pred_mean") |>
  dplyr::select(month, value, geometry) |>
  split(~month) |>
  lapply(function(x) {
    x |>
      dplyr::select(value, geometry) |>
      st_rasterize() |>
      rast() |>
      mask(sim[[1]])
  })

length(a)
plot(a[[2]])

unique(a$month)

i <- st_read("~/OneDrive - University of Glasgow/PhD/0_simulations/post_samples/anopheles/i_50_an_df.geojson") |>
  dplyr::filter(variable == "pred_mean") |>
  dplyr::select(month, value, geometry) |>
  split(~month) |>
  lapply(function(x) {
    x |>
      dplyr::select(value, geometry) |>
      st_rasterize() |>
      rast() |>
      mask(sim[[1]])
  })

sim_values <- sim[[1]][!is.na(sim[[1]])]
i_values <- i[!is.na(i)]
a_values <- a[!is.na(a)]

n_breaks <- 6
breaks_sim <- quantile(sim_values, probs = seq(0, 1, length.out = n_breaks), na.rm = TRUE)
breaks_a <- quantile(a_values, probs = seq(0, 1, length.out = n_breaks), na.rm = TRUE)
breaks_i <- quantile(i_values, probs = seq(0, 1, length.out = n_breaks), na.rm = TRUE)

sim_matrix <- matrix(c(breaks_sim[-length(breaks_sim)],
                           breaks_sim[-1],
                           seq_along(breaks_sim[-1]) -1),
                         ncol = 3, byrow = FALSE)

sim_classified <- classify(sim[[1]], sim_matrix)

a_matrix <- matrix(c(breaks_a[-length(breaks_a)],
                     breaks_a[-1],
                     seq_along(breaks_a[-1]) -1),
                   ncol = 3, byrow = FALSE)
a_matrix[1, 1] <- 0

a_classified <- classify(a,a_matrix)

i_matrix <- matrix(c(breaks_i[-length(breaks_i)],
                     breaks_i[-1],
                     seq_along(breaks_i[-1]) -1),
                   ncol = 3, byrow = FALSE)

i_matrix[1, 1] <- 0

i_classified <- classify(i,i_matrix)

n_classes <- length((unique(values(sim_classified, na.rm = TRUE))))
w_size <- 15
w_matrix <- matrix(1, nrow= w_size, ncol= w_size)

calculate_frequency <- function(classified_data) {
  lapply(0:(n_classes - 1), function(class) {
    dummy_class <- (classified_data == class) * 1
    freq <- focal(dummy_class, w=w_matrix, fun=sum, na.rm=TRUE) / (w_size^2)
    return(freq)
  })
}

sim_freq_list <- calculate_frequency(sim_classified)
a_freq_list <- calculate_frequency(a_classified)
i_freq_list <- calculate_frequency(i_classified)

calculate_squared_diff <- function(freq_list1, freq_list2) {
  lapply(seq_along(freq_list1), function(idx) {
    (freq_list1[[idx]] - freq_list2[[idx]])^2
  })
}

diff_a <- calculate_squared_diff(sim_freq_list, a_freq_list)
diff_i <- calculate_squared_diff(sim_freq_list, i_freq_list)

length(diff_a)
class(diff_a[[1]])

plot(diff_a[[1]])

sum_a <- Reduce("+", diff_a)
sum_i <- Reduce("+", diff_i)

plot(sum_a)
plot(sum_i)

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
