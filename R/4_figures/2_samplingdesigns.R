packs <- c("sf", "terra", "tidyverse", "raster", "rnaturalearth", "rnaturalearthdata", "tidyterra", "patchwork", "cowplot", "tmap")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> 
  st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
ss_points <- st_read('~/Documents/GitHub/OptimSS/data/3_sampling/sampling_designs.geojson') |>
  dplyr::filter(iteration == 7, sample_size == 15) |>
  mutate(
    gradient1 = case_when(
      scenario == "a" ~ 1,
      scenario == "b" ~ 2,
      scenario == "c" ~ 3,
      scenario == "d" ~ 1,
      scenario == "e" ~ 2,
      scenario == "f" ~ 3,
      scenario == "g" ~ 1,
      scenario == "h" ~ 2,
      scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      scenario == "a" ~ 1,
      scenario == "b" ~ 1,
      scenario == "c" ~ 1,
      scenario == "d" ~ 2,
      scenario == "e" ~ 2,
      scenario == "f" ~ 2,
      scenario == "g" ~ 3,
      scenario == "h" ~ 3,
      scenario == "i" ~ 3,
      TRUE ~ NA_real_
    )
  )

(ss_plot_effort <- ggplot() +
    geom_spatraster(data = as.factor(landcover)) +
    scale_fill_manual(
      values = c('#f6e8c3', '#80cdc1', '#018571', '#d8b365', 'gray'),
      labels = c("0" = "Oil palm", "1" = "Secondary forest", "2" = "Primary forest", "3" = "Other plantations", "4" = "Built-up"),
      na.value = NA,
      na.translate = FALSE,
      guide = guide_legend(nrow = 1)
    ) +
    geom_sf(data = boundary, fill = NA, color = "black") +
    geom_sf(data = ss_points, pch = 21, fill = 'red', alpha = .1, size = 2)+
  facet_grid(gradient2~gradient1, labeller = labeller(
    gradient1 = c("1" = "Lattice", "2" = "Stratified", "3" = "Random"),
    gradient2 = c("1" = "Static", "2" = "Rotational", "3" = "Dynamic")
  ))+
  theme_void(base_family = "Times New Roman", base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'bottom',
    text = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman"),
    plot.caption = element_text(family = "Times New Roman"),
    plot.margin = unit(rep(5, 4), "mm")
  ) +
  labs(title = NULL, 
       subtitle = NULL, 
       fill = "Land Use/Land Cover", 
       caption = "Columns: Spatial structure \nRows: Temporal flexibility")
)

ggsave('~/Documents/GitHub/OptimSS/data/5_figures/2_samplingdesigns.jpg', ss_plot_effort, dpi = 300, width = 8, height = 6)
