# Estimate parameters for the simulation
# Usage: Rscript R/1_simulation/1_EstimateParams.R [species]
# Example: Rscript R/1_simulation/1_EstimateParams.R anopheles
# This script estimates the parameters for the simulation of mosquito abundance 
# based on the provided species. It fits a spatial model using INLA and inlabru, 
# and saves the estimated parameters to a CSV file for further use in the simulation process.

### Load and install required packages
packs <- c("tidyverse", "INLA", "inlabru", "sf")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
### Load data
species <- commandArgs(trailingOnly = TRUE)
mos <- st_read(paste0('~/Documents/GitHub/OptimSS/data/1_raw/', species, "_raw.geojson"))
roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')

### Hyperparameters
if (species == "anopheles") {
  range <- 3000
  rho <- 0.935
} else if (species == "aedes") {
  range <- 1000
  rho <- 0.9
} 
### model fit
mesh <- fm_mesh_2d(roi, max.edge = c(2500, 5000), cutoff = 1000)
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, 0.01), prior.range = c(range, 0.1))
model <- bite_rate ~ Intercept(1) + 
  land_cover(Class, model = 'factor_contrast') + 
  field(geometry, model = matern, group = mon, 
        control.group = list(model = 'ar1', 
                             hyper = list(theta = list(prior = 'normal', param = c(rho, 0.001)))))

fit <- bru(model, mos, family = "nbinomial",
           options = list(control.family = list(link = "log"),
             control.inla = list(int.strategy = "eb")
           ))
# 
lc_summary <- fit$summary.random$land_cover[, c("ID", "mean", "sd")]
lc_summary$effect_type <- "land_cover"

intercept_summary <- data.frame(
  ID = "A Primary",
  mean = fit$summary.fixed["Intercept", "mean"],
  sd = fit$summary.fixed["Intercept", "sd"],
  effect_type = "land_cover"
)

combined_df <- rbind(intercept_summary, lc_summary)
write_csv(combined_df, paste0("~/Documents/GitHub/OptimSS/data/1_raw/", species, "_parameters.csv"))
