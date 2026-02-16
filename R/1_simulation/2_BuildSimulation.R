### Load and install required packages
packs <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("~/Documents/GitHub/OptimSS/R/1_simulation/helpers/mosquito_class.R")

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


compute_mean_surface <- function(boundary, lc_coefficients, spatial_range, rho, phi, land_cover, n_iterations = 100) {
    surfaces <- list()
    for (i in 1:n_iterations) {
        mosquito <- Mosquito(boundary = boundary,
                            lc_coefficients = lc_coefficients,
                            spatial_range = spatial_range,
                            rho = rho,
                            phi = phi,
                            land_cover = land_cover,
                            seed = i)
        surfaces[[i]] <- mosquito@simulated_surface
    }
    mean_surface <- Reduce("+", surfaces) / n_iterations
    return(list(surfaces = surfaces, mean_surface = mean_surface))
}

mean_surface <- compute_mean_surface(boundary = roi,
                                    lc_coefficients = beta,
                                    spatial_range = range,
                                    rho = rho,
                                    phi = phi,
                                    land_cover = landcover)


plot(mean_surface$mean_surface)

plot(mean_surface$surfaces[[100]])

any(sapply(1:(length(mean_surface$surfaces)-1), function(i) 
    identical(mean_surface$surfaces[[i]], mean_surface$surfaces[[i+1]])))

class(mean_surface)

# writeRaster(anopheles@simulated_surface, paste0(wd, '/data/20250331_sim_raster001.tif'), overwrite=TRUE)
writeRaster(mean_surface$mean_surface, '~/Documents/GitHub/OptimSS/data/1_raw/20260214_sim_raster001_mean.tif', overwrite=TRUE)
