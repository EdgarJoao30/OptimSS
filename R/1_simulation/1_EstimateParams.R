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
# mos <- st_read("~/Documents/GitHub/OptimSS/data/1_raw/aedes_raw.geojson")
roi <- st_read("~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson")

### Hyperparameters
if (species == "anopheles") {
  range <- 3000
  rho <- 0.9
} else if (species == "aedes") {
  range <- 1000
  rho <- 0.7
} 
### model fit
mesh <- fm_mesh_2d(roi, max.edge = c(2500, 5000), cutoff = 1000)

# prior.sigma = c(0.5, 0.01) says there is a 1% chance the marginal SD is above 0.5.
# prior.range = c(3000, 0.1) says there is a 10% chance the range is below 3000 m.
# if Prange or Psigma is NA, then range and sigma are used as fixed values and are not estimated.
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, NA), prior.range = c(range, NA))
# theta is a transformation of rho, so we set the prior on 
# theta to reflect our prior belief about rho. theta <- log((1 + rho) / (1 - rho))
theta <- log((1 + rho) / (1 - rho))
model <- bite_rate ~ Intercept(1) + 
  land_cover(Class, model = 'factor_contrast') + 
  field(geometry, model = matern, group = mon, 
        control.group = list(model = 'ar1', 
                             # param = c(theta, 10) says there is a 10% chance the autocorrelation is above rho.
                             # here 10 is the precision, which is the inverse of the variance. 
                             # A higher precision means a stronger belief in the prior, 
                             # while a lower precision means a weaker belief.
                             hyper = list(theta = list(prior = 'normal', param = c(theta, 10), initial = theta, fixed = TRUE))))

fit <- bru(model, mos, family = "nbinomial",
           control.family = list(
             link = "log", 
             variant = 0,
             hyper = list(size = list(prior = "normal", param = c(10, 1), initial = 10, fixed = TRUE))
           ),
           options = list(
             control.inla = list(int.strategy = "eb")
           ))

# Extract summary 
lc_summary <- fit$summary.random$land_cover[, c("ID", "mean", "sd")]
lc_summary$effect_type <- "land_cover"

intercept_summary <- data.frame(
  ID = "A Primary",
  mean = fit$summary.fixed["Intercept", "mean"],
  sd = fit$summary.fixed["Intercept", "sd"],
  effect_type = "intercept"
)

combined_df <- rbind(intercept_summary, lc_summary)
combined_df$exp <- exp(combined_df$mean)

combined_df <- combined_df |>
  mutate(
    mean_abundance = ifelse(
      effect_type == "land_cover",
      exp * filter(combined_df, ID == "A Primary")$exp,
      exp
    )
  )

write_csv(combined_df, paste0("~/Documents/GitHub/OptimSS/data/1_raw/", species, "_parameters.csv"))
