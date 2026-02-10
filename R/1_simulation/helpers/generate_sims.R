### Load and install required packages
packs <- c("tidyverse", "INLA", "inlabru", "sf", "fmesher", "terra", "stars", "raster")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)


generate_sims <- function(boundary = boundary,
                         landcover = landcover,
                         alpha = 2, 
                         theta = theta, 
                         rho = rho, 
                         beta = beta, 
                         k = k, # number of samples
                         sd_mu = sd_mu,
                         phi = phi,
                         seed = 1234) {
  
  points <- as.data.frame(landcover, xy = TRUE)
  # colnames(points)[3] <- 'class'
  # points$class <- factor(points$class)
  # points$class <- relevel(factor(points$class), ref = 2)
  points$class <- factor(points$Landcover_AllClass)
  points$class <- plyr::revalue(points$class, c("0"="C_Oil", "1"="B_Secondary", "2"="A_Primary", "3"="D_Plantation", "4"="E_Built"))
  points$class <- factor(points$class, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
  points <- st_as_sf(points, coords = c('x', 'y'), crs = 32650)
  mask <- rasterize(vect(boundary), landcover)
  # Mesh and true surface, units = meters
  mesh <- fm_mesh_2d(points, max.edge = c(2500, 5000), cutoff = 1000)
  spde <- inla.spde2.matern(mesh, alpha = alpha)
  Q <- inla.spde2.precision(spde, theta = theta)
  true_field <- inla.qsample(k, Q, seed = seed)
  points$field <- fm_evaluate(mesh, loc = points, field = true_field)
  # Compute AR1
  points$field_AR1 <- points$field
  for (j in 2:k) {
    points$field_AR1[, j] <- rho * points$field_AR1[, j - 1] + sqrt(1 - rho^2) * points$field[, j]
  }
  # Add regression covariates
  ccov <- factor(replicate(k, points$class))
  n <- nrow(points)
  
  
  beta0 <- beta[1]
  beta[1] <- 0
  #head(beta0 + beta[unclass(ccov)])
  
  mu <- beta0 + beta[unclass(ccov)] + points$field_AR1 + rnorm(n * k, 0, sd_mu)
  points$mu <- exp(mu)
  # Draw from nbinomial distribution
  generate_nbinomial <- function(x) {
    rnbinom(mu = x, n = 1, size = phi)
  }
  set.seed(seed)
  nbinomial_sample <- apply(points$mu, c(1, 2), generate_nbinomial)
  points$mosq <- nbinomial_sample
  points <- cbind(points, as.data.frame(points$mosq)) %>% dplyr::select(sample.1:sample.12)
  # Convert points to raster
  surface <- st_rasterize(points) %>% rast() %>% terra::mask(mask)
  names(surface) <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
  
  return(surface)
}



