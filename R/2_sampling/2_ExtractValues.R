# This script extracts the simulated mosquito abundance values and corresponding land cover classes at the sampling locations defined in the sampling design. 
# It saves the results as a GeoJSON file for each species.
# Usage: Rscript R/2_sampling/2_ExtractValues.R [species]
# Example: Rscript R/2_sampling/2_ExtractValues.R anopheles
### Load and install required packages
packs <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("./helpers/extract_values.R")

roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
sampling <- st_read('~/Documents/GitHub/OptimSS/data/3_sampling/sampling_designs.geojson')
species <- commandArgs(trailingOnly = TRUE)
simulation <- rast(paste0('~/Documents/GitHub/OptimSS/data/2_simulation/', species, '_sim.tif'))

values <- extract_values(sampling, landcover, simulation)

st_write(values, paste0("~/Documents/GitHub/OptimSS/data/3_sampling/", species, "_sampling.geojson"), append = F)
