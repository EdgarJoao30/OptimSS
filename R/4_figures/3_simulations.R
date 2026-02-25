packs <- c("sf", "terra", "tidyverse", "tidyterra", "patchwork", "cowplot", "extrafont")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
loadfonts(device = "pdf")
boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> 
  st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
anopheles <- rast('~/Documents/GitHub/OptimSS/data/2_simulation/anopheles_sim.tif')
aedes <- rast('~/Documents/GitHub/OptimSS/data/2_simulation/aedes_sim.tif')
anopheles_params <- read.csv('~/Documents/GitHub/OptimSS/data/1_raw/anopheles_parameters.csv')
aedes_params <- read.csv('~/Documents/GitHub/OptimSS/data/1_raw/aedes_parameters.csv')

(anoph_plot <- ggplot() +
  geom_spatraster(data = anopheles[[c('jan', 'jun', 'dec')]] + 1) +
  geom_sf(data = boundary, fill = NA, color = "black") +
  facet_wrap(
    ~lyr, ncol = 1,
    labeller = as_labeller(
      c(jan = "January", jun = "June", dec = "December")
    )
  ) +
  scale_fill_gradient2(
    low = "white", 
    high = "red",
    na.value = NA,
    trans = "log",
    breaks = c(1, 5, 50),
    labels = scales::label_number(accuracy = 1)
  ) +
  theme_void(base_family = 'Times New Roman', base_size = 12) +
  labs(title = expression(italic("(a) Anopheles balabacensis") ~ "simulation"), fill = 'Abundance',
       caption = paste0('Mean abundance: \
       Primary forest: ', round(anopheles_params$mean_abundance[anopheles_params$ID == "A Primary"], 2), ', ', 
       'Secondary forest: ', round(anopheles_params$mean_abundance[anopheles_params$ID == "B Secondary"], 2), ', ', '\n', 
       'Oil Palm: ', round(anopheles_params$mean_abundance[anopheles_params$ID == "C Oil palm"], 2), ', ', 
       'Other plantations: ', round(anopheles_params$mean_abundance[anopheles_params$ID == "D Plantation"], 2), ', ', 
       'Built-up: ', round(anopheles_params$mean_abundance[anopheles_params$ID == "E Built-up"], 2), ', ', '\n', 
       'Temporal correlation: 0.9, \
       Spatial range: 3000 meters, \ 
       Overdispersion: 10')) +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5), # Center the title
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman"),
    plot.caption = element_text(family = "Times New Roman"),
    plot.margin = unit(rep(5, 4), "mm")
  ) )

(aedes_plot <- ggplot() +
    geom_spatraster(data = aedes[[c('jan', 'jun', 'dec')]] + 1) +
    geom_sf(data = boundary, fill = NA, color = "black") +
    facet_wrap(
      ~lyr, ncol = 1,
      labeller = as_labeller(
        c(jan = "January", jun = "June", dec = "December")
      )
    ) + 
    scale_fill_gradient2(
      low = "white", 
      high = "red",
      na.value = NA,
      trans = "log",
      breaks = c(1, 5, 50),
      labels = scales::label_number(accuracy = 1)
    ) +
    theme_void(base_family = 'Times New Roman', base_size = 12) +
    labs(title = expression(italic("(b) Aedes albopictus") ~ "simulation"), fill = 'Abundance',
         caption = paste0('Mean abundance: \
       Primary forest: ', round(aedes_params$mean_abundance[aedes_params$ID == "A Primary"], 2), ', ', 
       'Secondary forest: ', round(aedes_params$mean_abundance[aedes_params$ID == "B Secondary"], 2), ', ', '\n', 
       'Oil Palm: ', round(aedes_params$mean_abundance[aedes_params$ID == "C Oil palm"], 2), ', ', 
       'Other plantations: ', round(aedes_params$mean_abundance[aedes_params$ID == "D Plantation"], 2), ', ', 
       'Built-up: ', round(aedes_params$mean_abundance[aedes_params$ID == "E Built-up"], 2), ', ', '\n', 
       'Temporal correlation: 0.7, \
       Spatial range: 1000 meters, \ 
       Overdispersion: 10')) +
    theme(
      legend.position = 'bottom',
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5), # Center the title
      text = element_text(family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      legend.title = element_text(family = "Times New Roman"),
      strip.text = element_text(family = "Times New Roman"),
      plot.caption = element_text(family = "Times New Roman"),
      plot.margin = unit(rep(5, 4), "mm")
    ) )

(mosquito_plot <- anoph_plot | aedes_plot)

ggsave('~/Documents/GitHub/OptimSS/data/5_figures/3_simulations.jpg', mosquito_plot, dpi = 300, width = 7, height = 7)
 