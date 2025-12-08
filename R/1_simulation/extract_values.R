#################################################### Extract values from habitat and simulations ###############################################################
library(sf)
library(terra)

extract_values <- function(df, lc, sims) {
  df_h <- terra::extract(lc, df) 
  df$land_cover <- df_h[,2]
  
  dfs <- list()
  
  for (m in 1:12) {
    df_m <- df %>% dplyr::filter(month == m)
    r_m <- sims[[m]]
    
    df_e <- terra::extract(r_m, df_m)
    colnames(df_e) <- c('ID', 'sim_anoph')
    df_m$sim_anoph <- df_e$sim_anoph
    dfs[[m]] <- df_m
  }
  
  df_all <- do.call(rbind, dfs)
  
  return(df_all)
}




wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |>
  st_union() |>
  st_transform(crs = 32650)
# Land cover
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif'))
df <- st_read(paste0(wd, '/data/test/20250728_points_sampling_scenarios_nobuffer.geojson'))
sim <- rast(paste0(wd, '/data/20250331_sim_raster001.tif'))

values <- extract_values(df, landcover, sim)
# 
# st_write(values, paste0(wd, '/data/test/20250728_points_sampling_scenarios_alldata.geojson'), append = F)
