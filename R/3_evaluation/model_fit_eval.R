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

species <- "aedes"
i <- 1
s <- 5

if (species == "anopheles") {
  w_size <- 13
} else if (species == "aedes") {
  w_size <- 5
} 

print(paste0('Running species: ', species))
print(paste0('Running sample size: ', s))
print(paste0('Running iteration: ', i))
# Landcover
boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
# Simulations
simulations <- rast(paste0('~/Documents/GitHub/OptimSS/data/2_simulation/', species, '_sim.tif'))
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


results_list <- lapply(scenarios, function(x) {
  # x <- 'b'
  sample <- samples %>% dplyr::filter(scenario == x)
  # Drop unused levels from the factor
  sample$cat_500m <- droplevels(as.factor(sample$cat_500m))
  # Set the first class as the baseline - Should be most frequent class but INLA orders them alphabetically
  first_class <- names(table(sample$cat_500m))[1]
  sample$cat_500m <- relevel(as.factor(sample$cat_500m), ref = first_class)
  missing_classes <- setdiff(levels(simulations$cat_500m), levels(sample$cat_500m))
  # 1. Determine the number of unique land cover classes
  n_levels <- length(unique(na.omit(sample$cat_500m)))

  # 2. Build the formula dynamically
  if (n_levels > 1) {
    model <- sim ~ Intercept(1) + 
      land_cover(cat_500m, model = 'factor_contrast') +
      field(geometry, model = matern, group = month, 
            control.group = list(model = 'ar1', 
                                hyper = list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))))
  } else {
    # If only one level, exclude the land_cover term entirely
    model <- sim ~ Intercept(1) + 
      field(geometry, model = matern, group = month, 
            control.group = list(model = 'ar1', 
                                hyper = list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))))
  }

  fit <- bru(model, sample, family = "nbinomial",
             options = list(control.family = list(link = "log"), 
               control.inla = list(int.strategy = "eb")
             ))
   
  summ <- summarize_model(fit, scenario = x, n_levels = n_levels) |>
    mutate(species = species, iteration = i, sample_size = s) |> 
    dplyr::select(species, scenario, sample_size, iteration, effect, mean) |>
    pivot_wider(names_from = effect, values_from = mean) 

  summ <- summ |>
    mutate(across(matches("^[A-E]_"), list(abundance = ~exp(Intercept + .x)), .names = "{.col}_abundance")) |>
    mutate(Intercept_abundance = exp(Intercept))
  
  summ <- summ |>
    rename(!!first_class := Intercept, !!paste0(first_class, "_abundance") := Intercept_abundance)

  # Add missing classes with NA values
  for (missing_class in missing_classes) {
    summ[[missing_class]] <- NA
    summ[[paste0(missing_class, "_abundance")]] <- NA
  }
  
  summ <- summ |>
    dplyr::select(species, scenario, sample_size, iteration, 
                  A_Primary, B_Secondary, C_Oil, D_Plantation, E_Built,
                  A_Primary_abundance, B_Secondary_abundance, C_Oil_abundance, D_Plantation_abundance, E_Built_abundance,
                  `size for the nbinomial observations (1/overdispersion)`, `Range for field`, `Stdev for field`, `GroupRho for field`)

  pred <- predict(fit, sim_long,
    if (n_levels > 1) {
      ~ list(mu = exp(Intercept + land_cover + field))
    } else {
      ~ list(mu = exp(Intercept + field))
    },
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

  metrics <- compute_metrics(sim, mean_rasters, w_size = w_size) |>
    as.data.frame() |>
    mutate(species = species, scenario = x, iteration = i, sample_size = s)

  all_metrics <- cbind(summ, metrics  |> dplyr::select(-species, -scenario, -iteration, -sample_size))
  return(all_metrics)
})

metrics_df <- results_list |> bind_rows()
