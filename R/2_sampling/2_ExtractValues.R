### Load and install required packages
packs <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("./helpers/extract_values.R")

roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
simulation <- rast('~/Documents/GitHub/OptimSS/data/2_simulation/20250331_sim_raster001.tif')
sampling <- st_read('~/Documents/GitHub/OptimSS/data/3_sampling/20250429_points_sampling_scenarios_nobuffer.geojson')

test_values <- extract_values(sampling, landcover, simulation)

# st_write(test_values, paste0(wd, '/data/20250331_points_sampling_scenarios_alldata.geojson'), append = F)
