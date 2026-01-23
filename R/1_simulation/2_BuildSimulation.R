### Load and install required packages
packs <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("./helpers/mosquito_class.R")

roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')
# Land cover
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')

#change this to pass a csv file with these values
beta <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367) # Primary, Secondary, Oil Palm, Plantation, Built-up
range <- 3000
rho <- 0.9353888
phi <- 10

anopheles <- Mosquito(boundary = roi,
                      lc_coefficients = beta, 
                      spatial_range = range,
                      rho = rho,
                      phi = phi,
                      land_cover = landcover,
                      seed = 1234)

plot(anopheles@simulated_surface)

# writeRaster(anopheles@simulated_surface, paste0(wd, '/data/20250331_sim_raster001.tif'), overwrite=TRUE)
