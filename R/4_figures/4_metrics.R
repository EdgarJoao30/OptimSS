packs <- c("tidyverse", "RColorBrewer", "patchwork", "sf", "terra", "tidyterra", "landscapemetrics", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/raster_management.R')

# Samples
species <- "anopheles"
i <- 15
s <- 15
design <- "e"
query_sql <- sprintf(
  "SELECT * FROM %s_sampling WHERE iteration = %s AND sample_size = %s AND scenario = '%s'",
  species, i, s, design
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
sample <- samples %>% dplyr::filter(scenario == design)
# Landcover
boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
# Simulations
sim <- rast('~/Documents/GitHub/OptimSS/data/2_simulation/anopheles_sim.tif')
sim_classified <- classify_raster_by_quantiles(sim)
simulations <- sim
simulations <- st_as_sf(as.data.frame(simulations, xy = TRUE), crs = 32650, coords = c('x', 'y'))
lc_values <- terra::extract(landcover, simulations, xy = TRUE)
simulations$cat_500m <- round(lc_values$Landcover_AllClass)
simulations$cat_500m <- plyr::revalue(as.character(simulations$cat_500m), 
                                      c("0" = "C_Oil", "1" = "B_Secondary", "2" = "A_Primary", "3" = "D_Plantation", "4" = "E_Built"))
simulations$cat_500m <- factor(simulations$cat_500m, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
sim_long <- simulations %>%
  pivot_longer(cols = -c(cat_500m:geometry), names_to = 'month', values_to = 'sim') %>%
  mutate(month = match(month, tolower(month.abb)))

# Landscape metrics
w_size <- 13
w_matrix <- matrix(1, nrow = w_size, ncol = w_size)
sim_lm <- window_lsm(sim_classified, window = w_matrix, 
                     what = c("lsm_l_joinent", "lsm_l_contig_mn"))

# INLA mesh and model definition
mesh <- fm_mesh_2d(simulations, max.edge = c(2500, 5000), cutoff = 1000)
# range being 1/10 the diagonal distance of the study area in meters
range <- sqrt((st_bbox(boundary)[2] - st_bbox(boundary)[1])^2 + (st_bbox(boundary)[4] - st_bbox(boundary)[3])^2) / 10
# sigma being 1/10 the standard deviation of the response variable in the sample data
sigma <- sd(samples$sim, na.rm = TRUE) / 10
# prior.sigma = c(1.5, 0.01) says there is a 1% chance the marginal SD is above 1.5.
# prior.range = c(8000, 0.9) says there is a 90% chance the range is below 8000 m.
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(sigma, 0.01), prior.range = c(range, 0.9))
# model definition
model <- sim ~ Intercept(1) + 
      land_cover(cat_500m, model = 'factor_contrast') +
      field(geometry, model = matern, group = month, 
            control.group = list(model = 'ar1', 
                                hyper = list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))))
# fit model 
fit <- bru(model, sample, family = "nbinomial",
           options = list(control.family = list(link = "log"), 
             control.inla = list(int.strategy = "eb")
           ))
# predict
pred <- predict(fit, sim_long,
  ~ list(mu = exp(Intercept + land_cover + field)),
  n.samples = 1000
)
# to raster
pred_rasters <- lapply(1:12, function(m) {
  pred$mu |>
    dplyr::select(month, mean, q0.025, q0.975) |>
    dplyr::filter(month == m) |>
    st_as_sf(geometry = "geometry") |>
    st_rasterize() |>
    rast() |>
    mask(landcover)
})
pred_mean <- lapply(pred_rasters, function(r) r[['mean_mean']]) |> rast()
pred_classified <- classify_raster_by_quantiles(pred_mean)
pred_lm <- window_lsm(pred_classified, window = w_matrix, 
                      what = c("lsm_l_joinent", "lsm_l_contig_mn"))

# plot sim and sim classified side by side using ggplot and tidyterra

plot(sim_lm$layer_1$lsm_l_joinent |> mask(sim[[1]]))
plot(pred_lm$layer_1$lsm_l_joinent |> mask(sim[[1]]))

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

plot(c(pred_mean[[1]], 
       pred_classified[[1]], 
       pred_lm$layer_1$lsm_l_joinent |> mask(sim[[1]]),
       pred_lm$layer_1$lsm_l_contig_mn |> mask(sim[[1]])))
