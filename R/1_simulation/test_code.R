packages <- c("sf", "terra", "tidyverse", "ggridges", "tidyterra", "stars", "raster", "INLA", "inlabru", "fmesher")

# Check if packages are installed, and install if not
installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load the packages
library(sf)
library(terra)
library(tidyverse)
library(ggridges)
library(tidyterra)
library(stars)
library(raster)
library(INLA)
library(inlabru)
library(fmesher)
source("./R/mosquito_class.R")
source("./R/extract_values.R")

wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |>
  st_union() |>
  st_transform(crs = 32650)
# Land cover
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif'))
df <- st_read(paste0(wd, '/data/20250129_points_sampling_scenarios_nobuffer.geojson'))

beta <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367) # Primary, Secondary, Oil Palm, Plantation, Built-up
range <- 3000
rho <- 0.9353888
phi <- 10

anopheles <- Mosquito(lc_coefficients = beta, 
                      spatial_range = range,
                      rho = rho,
                      phi = phi,
                      land_cover = landcover,
                      seed = 1234)

plot(anopheles@simulated_surface)

# writeRaster(anopheles@simulated_surface, paste0(wd, '/data/20250331_sim_raster001.tif'), overwrite=TRUE)

test_values <- extract_values(df, anopheles@land_cover, anopheles@simulated_surface)

# st_write(test_values, paste0(wd, '/data/20250331_points_sampling_scenarios_alldata.geojson'), append = F)

lc <- as.data.frame(anopheles@land_cover, xy=T)
sims <- as.data.frame(anopheles@simulated_surface, xy=T)
sims$lc <- as.factor(lc[,3])
sims <- gather(sims, key = 'month', value = 'value', -c(x, y, lc)) %>% st_as_sf(coords = c('x', 'y'))

(ridge_plot <- ggplot(sims %>% dplyr::filter(month == 'mar'), aes(x = value, y = lc, fill = lc)) +
    geom_density_ridges(scale = 1, alpha = 0.7) +
    scale_fill_manual(values = c('darkgreen', '#9fa86a', 'green', '#def016', 'gray')) +
    labs(
      x = "Biting Rate",
      y = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
)

