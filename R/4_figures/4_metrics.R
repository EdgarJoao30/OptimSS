packs <- c("tidyverse", "RColorBrewer", "patchwork", "sf", "terra", "tidyterra")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

sim <- rast('~/Documents/GitHub/OptimSS/data/2_simulation/anopheles_sim.tif')
sim_classified <- classify_raster_by_quantiles(sim)
w_size = 13
w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
sim_lm <- window_lsm(sim_classified, window = w_matrix, 
                       what = c("lsm_l_joinent", "lsm_l_contig_mn"))

# plot sim and sim classified side by side using ggplot and tidyterra

plot(sim_lm$layer_1$lsm_l_joinent |> mask(sim[[1]]))

ggplot() +
  geom_spatraster(data = c(sim[[1]], sim_classified[[1]])) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_viridis_c(na.value = NA) +
  ggtitle("Original Simulation") +
  theme_void()

plot(c(sim[[1]], 
       sim_classified[[1]], 
       sim_lm$layer_1$lsm_l_joinent |> mask(sim[[1]]),
       sim_lm$layer_1$lsm_l_contig_mn |> mask(sim[[1]])))


# Samples
species <- "anopheles"
i <- 1
s <- 15
scenario <- "e"
query_sql <- sprintf(
  "SELECT * FROM %s_sampling WHERE iteration = %s AND sample_size = %s AND scenario = '%s'",
  species, i, s, scenario
)
samples <- st_read(
  dsn = paste0('~/Documents/GitHub/OptimSS/data/3_sampling/', species, '_sampling.geojson'),
  query = query_sql,
  quiet = TRUE # Optional: silences the driver output
)
samples <- samples |>
  mutate(
    cat_500m = factor(
      recode(as.character(land_cover), 
             "0" = "C_Oil", "1" = "B_Secondary", "2" = "A_Primary", 
             "3" = "D_Plantation", "4" = "E_Built"),
      levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built")
    ) 
  ) |> 
  rename(geometry = `_ogr_geometry_`)
