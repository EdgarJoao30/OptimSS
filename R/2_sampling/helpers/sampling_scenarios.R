#################################################### Generate sampling scenarios ###############################################################
# This function generates sampling scenarios based on the specified parameters. 
# It creates samples and grids for each combination of model, sample size, and iteration.
# The function takes the following arguments:
# - roi: A spatial object representing the region of interest.
# - iterations: The number of iterations to perform for each combination of model and sample size.
# - sample_sizes: A vector of sample sizes to use for generating samples.
# - models: A vector of model identifiers to specify the sampling strategy.
# - crs: The coordinate reference system to use for the spatial data.
packs <- c("sf", "tidyverse")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/2_sampling/helpers/sampling_functions.R')

generate_sampling_scenarios <- function(roi,
                    iterations = 100,
                    sample_sizes = c(5, 10, 15, 25, 50),
                    models = c('a','b','c','d','e','f','g','h','i'),
                    crs = 32650) {
  samples <- list()
  grids <- list()
  for (m in seq_along(models)) {
  model_s_temp <- list()
  model_g_temp <- list()
  for (n in seq_along(sample_sizes)) {
    samples_temp <- list()
    grids_temp <- list()
    for (i in seq_len(iterations)) {
    if (models[m] == 'a') {
      # Lattice - static sampling design
      temp <- uniform_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'b') {
      # Stratified - static sampling design
      temp <- stratified_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'c') {
      # Random - static sampling design
      temp <- random_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'd') {
      # Lattice - rotational sampling design
      ss1 <- round(sample_sizes[n] / 2)
      ss2 <- sample_sizes[n] - ss1
      temp1 <- uniform_sample(roi, sample_size = ss1, fixed = TRUE)
      temp2 <- uniform_sample(roi, sample_size = ss2, fixed = FALSE)
      s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = crs)
      s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = crs)
      s <- rbind(s1, s2)
      g1 <- temp1[[2]]
      g2 <- temp2[[2]]
      g <- rbind(g1, g2)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'e') {
      # Stratified - rotational sampling design
      ss1 <- round(sample_sizes[n] / 2)
      ss2 <- sample_sizes[n] - ss1
      temp1 <- stratified_sample(roi, sample_size = ss1, fixed = TRUE)
      temp2 <- stratified_sample(roi, sample_size = ss2, fixed = FALSE)
      s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = crs)
      s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = crs)
      s <- rbind(s1, s2)
      g1 <- temp1[[2]]
      g2 <- temp2[[2]]
      g <- rbind(g1, g2)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'f') {
      # Random - rotational sampling design
      ss1 <- round(sample_sizes[n] / 2)
      ss2 <- sample_sizes[n] - ss1
      temp1 <- random_sample(roi, sample_size = ss1, fixed = TRUE)
      temp2 <- random_sample(roi, sample_size = ss2, fixed = FALSE)
      s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = crs)
      s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = crs)
      s <- rbind(s1, s2)
      g1 <- temp1[[2]]
      g2 <- temp2[[2]]
      g <- rbind(g1, g2)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'g') {
      # Lattice - dynamic sampling design
      temp <- uniform_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'h') {
      # Stratified - dynamic sampling design
      temp <- stratified_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else if (models[m] == 'i') {
      # Random - dynamic sampling design
      temp <- random_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
      print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
    } else {
      next
    }
    if (!exists("s")) {
      s <- st_as_sf(temp[[1]], coords = c('x', 'y'), crs = crs)
      g <- temp[[2]]
    }
    s$iteration <- i
    g$iteration <- i
    s$sample_size <- sample_sizes[n]
    g$sample_size <- sample_sizes[n]
    s$scenario <- models[m]
    g$scenario <- models[m]
    samples_temp[[i]] <- s
    grids_temp[[i]] <- g
    rm(s)
    rm(g)
    } # end iterations
    model_s_temp[[n]] <- do.call(rbind, samples_temp)
    model_g_temp[[n]] <- do.call(rbind, grids_temp)
  } # end sample_sizes
  samples[[m]] <- do.call(rbind, model_s_temp)
  grids[[m]] <- do.call(rbind, model_g_temp)
  } # end designs
  samples <- do.call(rbind, samples)
  grids <- do.call(rbind, grids)
  grids <- grids[, c('id', 'cat', 'month', 'iteration', 'sample_size', 'scenario', 'geometry')]
  return(list(samples = samples, grids = grids))
}
