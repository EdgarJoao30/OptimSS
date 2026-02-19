#################################################### Model fitting ###############################################################
packs <- c("sf", "terra", "tidyverse", "raster", "INLA", "inlabru", "fmesher", "ModelMetrics", "parallel")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/3_evaluation/helpers/metrics.R')
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
boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/20240312_ROI_4326.shp') |> st_union() |> st_transform(crs = 32650) 
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
samples <- st_read(paste0('~/Documents/GitHub/OptimSS/data/3_sampling/', species, '_sampling.geojson'))
samples$cat_500m <- plyr::revalue(as.character(samples$land_cover), 
                                  c("0" = "C_Oil", "1" = "B_Secondary", "2" = "A_Primary", "3" = "D_Plantation", "4" = "E_Built"))
samples$cat_500m <- factor(samples$cat_500m, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
# Sub-sampling scenarios
scenarios <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
samples_list <- lapply(scenarios, function(x) {
  samples %>% dplyr::filter(iteration == i, scenario == x, sample_size == s)
})
names(samples_list) <- scenarios
# INLA mesh and model definition
mesh <- fm_mesh_2d(simulations, max.edge = c(2500, 5000), cutoff = 1000)
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, 0.01), prior.range = c(1000, 0.1))
# inlabru does not automatically remove the first level and absorb it into an intercept. 
# Instead, we can either use model="factor_full" without an intercept, or model="factor_contrast", 
# which does remove the first level.
# https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html
model <- sim ~ Intercept(1)  + 
  land_cover(cat_500m, model = 'factor_contrast') +
  field(geometry, model = matern, group = month, control.group = list(model = 'ar1'))

print(paste0('Model definition: ', model[3]))
start_time <- Sys.time()
fit_list <- lapply(samples_list, function(sample) {
  bru(model, sample, family = "nbinomial",
      options = list(control.family = list(link = "log"), 
                     #control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, waic = TRUE),
                     control.inla = list(int.strategy = "eb")
      ))
})
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat("Elapsed time:", elapsed_time, "seconds\n")

start_time <- Sys.time()
pred_list <- lapply(fit_list, function(fit) {
  predict(fit, sim_long,
          ~ {
            mu <- exp(Intercept + land_cover + field)
            
            list(
              mu = mu
            )
          },
          n.samples = 1000
  )
})
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat("Elapsed time:", elapsed_time, "seconds\n")

pred_raster <- lapply(1:9, function(scenarios) {
  lapply(1:12, function(m) {
    st_as_sf(pred_list[[scenarios]]$mu |> dplyr::filter(month == m), geometry = "geometry") |>
      st_rasterize() |>
      rast() |>
      mask(landcover)
  }) |> rast()
})

length(pred_raster)
plot(pred_raster[[9]]['mean_mean'])

class(pred_raster[[2]]['mean_mean'])
names(pred_raster[[11]])

start_time <- Sys.time()
metrics_list <- mclapply(seq_along(scenarios), function(idx) {
  metrics <- compute_metrics(sim, pred_raster[[idx]]['mean_mean'], w_size = 13)
  df <- as.data.frame(metrics)
  df$scenario <- scenarios[idx]
  df
}, mc.cores = 2)

metrics_df <- bind_rows(metrics_list)
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat("Elapsed time:", elapsed_time, "seconds\n")



# Extract mean and sd for fixed effects, random effects, and hyperparameters for each model in fit_list
results_df <- do.call(rbind, lapply(seq_along(fit_list), summarize_model))
results_df$iteration <- i
results_df$sample_size <- s

# Reset rownames for the final data frame
rownames(results_df) <- NULL

save(fit_list, file = paste0("/Volumes/Aedes/PhD/models/", "fitted_s", s, "_i", i ,".RData"))

print('DONE! with model fitting')

pred_list <- lapply(fit_list, function(fit) {
  predict(fit, sim_long,
          ~ {
            mu <- exp(Intercept + land_cover + field)
            predicted <- rnbinom(n = nrow(sim_long), size = fit$summary.hyperpar$mean[1], mu = mu)
            rmse <- ModelMetrics::rmse(sim, predicted)

            list(
              mu = mu,
              predicted = predicted,
              rmse = rmse
            )
          },
          n.samples = 1000
  )
})

# Add RMSE values to results_df

results_df <- lapply(seq_along(pred_list), function(idx) {
  means <- pred_list[[idx]]$rmse$mean
  sds <- pred_list[[idx]]$rmse$sd
  
  df <- results_df %>% dplyr::filter(Scenario == scenarios[idx])
  
  row <- df[1, ]
  row$Effect <- "rmse"
  row$Mean <- means
  row$SD <- sds
  row$Type <- "Error"
  
  rbind(df, row)
}) %>% bind_rows()

write.csv(results_df, paste0("/Volumes/Aedes/PhD/results/", 'summary_s', s, '_i', i, '.csv'))

print('DONE! predicting!')

samp_list <- lapply(fit_list, function(fit) {
  generate(fit, sim_long,
          ~ exp(Intercept + land_cover + field),
          n.samples = 1000
  )
})

save(samp_list, file = paste0("/Volumes/Aedes/PhD/post_samples/", "post_s", s, "_i", i ,".RData"))

print('DONE! sampling!')