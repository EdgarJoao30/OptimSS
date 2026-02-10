# This script processes the model results and simulation raster data to create a GeoJSON file containing 
# the predicted mean, standard deviation, and RMSE for Anopheles mosquitoes for a specific month. 
# The resulting GeoJSON file can be used for spatial analysis and visualization in GIS software.
packs <- c("tidyverse", "terra", "sf", "ModelMetrics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)


species <- commandArgs(trailingOnly = TRUE)[1]
sampling_strategy <- commandArgs(trailingOnly = TRUE)[2]
sampling_size <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations' 
ssdwd <- paste0("/Volumes/", R.utils::capitalize(species), "/PhD")
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |>
  st_union() |>
  st_transform(crs = 32650)

model_result <- paste0(ssdwd, "/post_samples/lite/post_", sampling_strategy, "_", sampling_size, ".RData")
sim_raster <- paste0(wd, '/data/20250331_sim_raster001.tif')

preprocess_df <- function(model_result, sim_raster, m) {
  load(model_result)
  mosquito <- sampled_data_list
  rm(sampled_data_list)
  sim <- rast(sim_raster) |>
    as.data.frame(xy = TRUE) |>
    st_as_sf(crs = 32650, coords = c('x', 'y')) |>
    pivot_longer(cols = -c(geometry), names_to = 'month', values_to = 'sim') |>
    mutate(month = match(month, tolower(month.abb))) 
  
  mosquito <- mosquito[[names(mosquito)]][sim$month == m, ]
  sim <- sim %>% filter(month == m)
  
  mosquito_se <- (sim$sim - mosquito)^2
  
  mosquito <- sim %>% 
    mutate(
      pred_mean = apply(mosquito, 1, mean, na.rm = TRUE),
      pred_sd = apply(mosquito, 1, sd, na.rm = TRUE),
      rmse = sqrt(apply(mosquito_se, 1, mean, na.rm = TRUE))
    )
  
  # mat <- mosquito
  # truth <- as.matrix(sim$sim)
  # truth_mat <- matrix(rep(truth, times = ncol(mat)), nrow = nrow(mat), ncol = ncol(mat))
  # rmse_vector <- sqrt(rowMeans((mat - truth_mat)^2))
  
  mosquito_long <- mosquito |>
    pivot_longer(cols = c(sim, pred_mean, pred_sd, rmse), names_to = "variable", values_to = "value")
    
  return(mosquito_long)
} 


df <- map_df(1:12, ~preprocess_df(model_result, sim_raster, .x))

st_write(df, paste0(wd, '/post_samples/', species, '/', sampling_strategy, '_', sampling_size, '_allmonths.geojson'))
