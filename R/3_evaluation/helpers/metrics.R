packs <- c("tidyverse", "terra", "sf", "stars", "landscapemetrics", "irr")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/local_metrics.R')
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/focal_metrics.R')
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/global_metrics.R')

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

compute_metrics <- function(sim, a, i, w_size = 13) {
    l <- local_metrics(sim, a)
    f <- focal_metrics(sim, a, w_size = w_size)
    g <- global_metrics(sim, a)
    list(local = l, focal = f, global = g)
}

# all_metrics <- compute_metrics(sim, a, i)
