#################################################### Extract values from habitat and simulations ###############################################################
# This script extracts the land cover and simulated values for each sampling point in the sampling scenarios. 
# The output is a geojson file with the sampling points and the extracted values for land cover and simulations.
### Load and install required packages
packs <- c("sf", "terra")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

extract_values <- function(df, lc, sims) {
  df_h <- terra::extract(lc, df) 
  df$land_cover <- df_h[,2]
  
  dfs <- list()
  
  for (m in 1:12) {
    df_m <- df %>% dplyr::filter(month == m)
    r_m <- sims[[m]]
    
    df_e <- terra::extract(r_m, df_m)
    colnames(df_e) <- c('ID', 'sim')
    df_m$sim <- df_e$sim
    dfs[[m]] <- df_m
  }
  
  df_all <- do.call(rbind, dfs)
  
  return(df_all)
}