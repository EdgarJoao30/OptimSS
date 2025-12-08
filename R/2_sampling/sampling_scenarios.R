#################################################### Generate sampling scenarios ###############################################################
library(sf)
library(tidyverse)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
source('~/Documents/GitHub/sampling_sims/R/sampling_functions.R')
roi <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_transform(crs = 32650) #%>% st_buffer(-10)

############## 
############## 
############## 
# iterations = 100
iterations = 1
# sample_sizes = c(5, 10, 15, 25, 50)
sample_sizes = c(300, 1000)
models = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')

samples <- list()
grids <- list()
for (m in 1:length(models)) {
  model_s_temp <- list()
  model_g_temp <- list()
  for (n in 1:length(sample_sizes)) {
    samples_temp <- list()
    grids_temp <- list()
    for (i in 1:iterations) {
      if (m ==1 ) {
        temp <- uniform_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 2) {
        temp <- stratified_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 3) {
        temp <- random_sample(roi, sample_size = sample_sizes[n], fixed = TRUE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 4) {
        ss1 <- round(sample_sizes[n] / 2)
        ss2 <- sample_sizes[n] - ss1
        temp1 <- random_sample(roi, sample_size = ss1, fixed = TRUE)
        temp2 <- uniform_sample(roi, sample_size = ss2, fixed = FALSE)
        
        s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = 32650)
        s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = 32650)
        s <- rbind(s1, s2)
        
        g1 <- temp1[[2]]
        g2 <- temp2[[2]]
        g <- rbind(g1, g2)
        
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 5) {
        ss1 <- round(sample_sizes[n] / 2)
        ss2 <- sample_sizes[n] - ss1
        temp1 <- random_sample(roi, sample_size = ss1, fixed = TRUE)
        temp2 <- stratified_sample(roi, sample_size = ss2, fixed = FALSE)
        
        s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = 32650)
        s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = 32650)
        s <- rbind(s1, s2)
        
        g1 <- temp1[[2]]
        g2 <- temp2[[2]]
        g <- rbind(g1, g2)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 6) {
        ss1 <- round(sample_sizes[n] / 2)
        ss2 <- sample_sizes[n] - ss1
        temp1 <- random_sample(roi, sample_size = ss1, fixed = TRUE)
        temp2 <- random_sample(roi, sample_size = ss2, fixed = FALSE)
        
        s1 <- st_as_sf(temp1[[1]], coords = c('x', 'y'), crs = 32650)
        s2 <- st_as_sf(temp2[[1]], coords = c('x', 'y'), crs = 32650)
        s <- rbind(s1, s2)
        
        g1 <- temp1[[2]]
        g2 <- temp2[[2]]
        g <- rbind(g1, g2)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 7) {
        temp <- uniform_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 8) {
        temp <- stratified_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else if (m == 9) {
        temp <- random_sample(roi, sample_size = sample_sizes[n], fixed = FALSE)
        print(paste0('Model: ', models[m], ' - ', 'Sample size: ', sample_sizes[n], ' - ', 'Iteration: ', i))
      } else {
        next  # Skip the iteration if the condition is not met
      }
      
      if(!exists("s")){
        s <- st_as_sf(temp[[1]], coords = c('x', 'y'), crs = 32650)
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
    }
    
    model_s_temp[[n]] <- do.call(rbind, samples_temp)
    model_g_temp[[n]] <- do.call(rbind, grids_temp)
    
  }
  
  samples[[m]] <- do.call(rbind, model_s_temp)
  grids[[m]] <- do.call(rbind, model_g_temp)
  
}
samples <- do.call(rbind, samples)
grids <- do.call(rbind, grids)
grids <- grids[, c('id', 'cat', 'month', 'iteration', 'sample_size', 'scenario', 'geometry')]

st_write(samples, paste0(wd, '/data/test/20250728_points_sampling_scenarios_nobuffer.geojson'), append = F)
st_write(grids, paste0(wd, '/data/test/20250728_grids_sampling_scenarios_nobuffer.geojson'), append = F)
