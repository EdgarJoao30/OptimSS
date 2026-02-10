### Load and install required packages
packs <- c("sf", "terra", "stars", "raster", "INLA", "inlabru", "fmesher")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source("~/Documents/GitHub/OptimSS/R/1_simulation/helpers/generate_sims.R")

# Define mosquito class
setClass(
  "Mosquito",
  slots = list(
    boundary = 'sf',
    lc_coefficients = "numeric",
    spatial_range = "numeric",
    rho = 'numeric',
    phi = 'numeric',
    theta = 'numeric',
    land_cover = 'SpatRaster',
    simulated_surface = 'SpatRaster',
    seed = 'numeric'
  )
)

setGeneric("_get_theta", function(object) {
  standardGeneric("_get_theta")
})

setMethod("_get_theta", "Mosquito", function(object) {
  alpha <- 2
  sigma <- 0.5
  variance <- sigma^2
  range <-object@spatial_range
  kappa <- sqrt(8 * (alpha - 1)) / range
  object@theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  return(object)  # Return the modified object 
})

setGeneric("_get_sims", function(object) {
  standardGeneric("_get_sims")
})

setMethod("_get_sims", "Mosquito", function(object) {
  result <- generate_sims(object@boundary,
                          object@land_cover,
                          theta = object@theta, 
                          rho = object@rho, 
                          beta = object@lc_coefficients, 
                          k = 12, 
                          sd_mu = 0.001,
                          phi = object@phi,
                          seed = 1234) 
  
  object@simulated_surface <- result
  return(object)  # Return the modified object 
})


# Constructor function for the Mosquito class
Mosquito <- function(boundary, lc_coefficients, spatial_range, rho, phi, land_cover, seed) {
  obj <- new("Mosquito", 
             boundary = boundary,
             lc_coefficients = lc_coefficients, 
             spatial_range = spatial_range,
             rho = rho,
             phi = phi,
             land_cover = land_cover,
             seed = seed)
  
  obj <- `_get_theta`(obj)  # Automatically apply the _get_theta method
  obj <- `_get_sims`(obj) 
  return(obj)
}


#beta <- c(0, 0.4567584, 1.5648494, 1.0986123, 1.8282377)
# beta <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367) # Primary, Secondary, Oil Palm, Plantation, Built-up
# range <- 3000
# rho <- 0.9353888
# phi <- 10
# 
# anopheles <- Mosquito(lc_coefficients = beta, 
#                       spatial_range = range,
#                       rho = rho,
#                       phi = phi,
#                       land_cover = landcover,
#                       seed = 1234)
# 
# plot(anopheles@simulated_surface)
# 
# # writeRaster(anopheles@simulated_surface, paste0(wd, '/data/20250331_sim_raster001.tif'), overwrite=TRUE)
# 
# test_values <- extract_values(df, anopheles@land_cover, anopheles@simulated_surface)
# 
# # st_write(test_values, paste0(wd, '/data/20250331_points_sampling_scenarios_alldata.geojson'), append = F)

