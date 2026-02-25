#################################################### Model fitting ###############################################################
packs <- c("sf", "terra", "tidyverse", "raster", "INLA", "inlabru", "fmesher", "ModelMetrics", "parallel")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/metrics.R')
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/summarize_model.R')
# user inputs
species <- commandArgs(trailingOnly = TRUE)[1]
i <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
s <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

species <- "anopheles"
i <- 1
s <- 15

print(paste0('Running species: ', species))
print(paste0('Running sample size: ', s))
print(paste0('Running iteration: ', i))
# Landcover
boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif', quiet = TRUE)
# Simulations
simulations <- rast(paste0('~/Documents/GitHub/OptimSS/data/2_simulation/', species, '_sim.tif'), quiet = TRUE)
sim <- simulations
simulations <- st_as_sf(as.data.frame(simulations, xy = TRUE), crs = 32650, coords = c('x', 'y'))
lc_values <- terra::extract(landcover, simulations, xy = TRUE)
simulations$cat_500m <- round(lc_values$Landcover_AllClass)
simulations$cat_500m <- plyr::revalue(as.character(simulations$cat_500m), 
                                      c("0" = "C_Oil", "1" = "B_Secondary", "2" = "A_Primary", "3" = "D_Plantation", "4" = "E_Built"))
simulations$cat_500m <- factor(simulations$cat_500m, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
sim_long <- simulations %>%
  pivot_longer(cols = -c(cat_500m:geometry), names_to = 'month', values_to = 'sim') %>%
  mutate(month = match(month, tolower(month.abb)))
# Samples
query_sql <- sprintf(
  "SELECT * FROM %s_sampling WHERE iteration = %s AND sample_size = %s",
  species, i, s
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
# Sub-sampling scenarios
scenarios <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
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
# inlabru does not automatically remove the first level and absorb it into an intercept. 
# Instead, we can either use model="factor_full" without an intercept, or model="factor_contrast", 
# which does remove the first level.
# https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html
model <- sim ~ Intercept(1)  + 
  land_cover(cat_500m, model = 'factor_contrast') +
  field(geometry, model = matern, group = month, 
    control.group = list(model = 'ar1',
                         # Sets rho = 0 with 90% chance of being positive
                         hyper = list(rho = list(prior = 'pc.cor1', param = c(0, 0.9))))
  )

results_list <- lapply(scenarios, function(x) {
  # x <- 'a'
  sample <- samples %>% dplyr::filter(scenario == x)
  # Drop unused levels from the factor
  sample$cat_500m <- droplevels(as.factor(sample$cat_500m))

  fit <- bru(model, sample, family = "nbinomial",
             options = list(control.family = list(link = "log"), 
               control.inla = list(int.strategy = "eb")
             ))
  
  summ <- summarize_model(fit, scenario = x) |>
    mutate(species = species, iteration = i, sample_size = s) |> 
    select(species, scenario, sample_size, iteration, effect, mean) |>
    pivot_wider(names_from = effect, values_from = mean) 

  pred <- predict(fit, sim_long,
    ~ list(mu = exp(Intercept + land_cover + field)),
    n.samples = 1000
  )

  pred_rasters <- lapply(1:12, function(m) {
    pred$mu |>
      dplyr::select(month, mean, q0.025, q0.975) |>
      dplyr::filter(month == m) |>
      st_as_sf(geometry = "geometry") |>
      st_rasterize() |>
      rast() |>
      mask(landcover)
  })

  mean_rasters <- lapply(pred_rasters, function(r) r[['mean_mean']]) |> rast()

  metrics <- compute_metrics(sim, mean_rasters, w_size = 13) |>
    as.data.frame() |>
    mutate(species = species, scenario = x, iteration = i, sample_size = s)

  all_metrics <- cbind(summ, metrics  |> dplyr::select(-species, -scenario, -iteration, -sample_size))
  return(all_metrics)
})

metrics_df <- results_list |> bind_rows()
