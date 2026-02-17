# This script simulates mosquito abundance surfaces based on the parameters estimated in 1_EstimateParams.R.
# It uses the Mosquito class defined in helpers/mosquito_class.R to generate simulated surfaces for a given species, 
# and saves the results as GeoTIFF files.
# Usage: Rscript R/1_simulation/2_BuildSimulation.R [species]
# Example: Rscript R/1_simulation/2_BuildSimulation.R anopheles
### Load and install required packages
packs <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("~/Documents/GitHub/OptimSS/R/1_simulation/helpers/mosquito_class.R")

roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')

species <- commandArgs(trailingOnly = TRUE)
# Order of land cover classes in files: Primary, Secondary, Oil Palm, Plantation, Built-up
beta <- read.csv(paste0('~/Documents/GitHub/OptimSS/data/1_raw/', species, '_parameters.csv'))$mean 

### Hyperparameters
if (species == "anopheles") {
  range <- 3000
  rho <- 0.9
} else if (species == "aedes") {
  range <- 1000
  rho <- 0.7
} 
phi <- 10 # size parameter for negative binomial distribution

mos_sim <- Mosquito(boundary = roi,
                      lc_coefficients = beta, 
                      spatial_range = range,
                      rho = rho,
                      phi = phi,
                      land_cover = landcover,
                      seed = 1234)

mean_sim <- compute_mean_surface(boundary = roi,
                                    lc_coefficients = beta,
                                    spatial_range = range,
                                    rho = rho,
                                    phi = phi,
                                    land_cover = landcover)

writeRaster(mos_sim@simulated_surface, paste0("~/Documents/GitHub/OptimSS/data/2_simulation/", species, '_sim.tif'), overwrite=TRUE)
writeRaster(mean_sim$mean_surface, paste0("~/Documents/GitHub/OptimSS/data/2_simulation/", species, '_sim__mean.tif'), overwrite=TRUE)
