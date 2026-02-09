library(terra)
library(sf)
library(stars)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(landscapemetrics)

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

lsm_l_pd(sim_classified)
freq(sim_classified)

w_matrix <- matrix(1, nrow=15, ncol=15)

sim_0 <- sim_classified
sim_0[sim_0 == 0] <- 11
sim_0[sim_0 < 11] <- 0
sim_0[sim_0 == 11] <- 1

sim_1 <- sim_classified
sim_1[sim_1 != 1] <- 0

sim_2 <- sim_classified
sim_2[sim_2 != 2] <- 0
sim_2[sim_2 == 2] <- 1

sim_3 <- sim_classified
sim_3[sim_3 != 3] <- 0
sim_3[sim_3 == 3] <- 1

sim_4 <- sim_classified
sim_4[sim_4 != 4] <- 0
sim_4[sim_4 == 4] <- 1

freq_window_0 <- focal(sim_0, w=w_matrix, fun=sum, na.rm=TRUE) / 225
freq_window_1 <- focal(sim_1, w=w_matrix, fun=sum, na.rm=TRUE) / 225
freq_window_2 <- focal(sim_2, w=w_matrix, fun=sum, na.rm=TRUE) / 225
freq_window_3 <- focal(sim_3, w=w_matrix, fun=sum, na.rm=TRUE) / 225
freq_window_4 <- focal(sim_4, w=w_matrix, fun=sum, na.rm=TRUE) / 225

plot(freq_window_0)
plot(freq_window_1)
plot(freq_window_2)
plot(freq_window_3)
plot(freq_window_4)
