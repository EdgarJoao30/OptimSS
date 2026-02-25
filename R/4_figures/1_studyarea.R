packs <- c("sf", "terra", "tidyverse", "raster", "rnaturalearth", "tidyterra", "patchwork", "cowplot", "tmap")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
data("World") ## get built-in data from the tmap package

boundary <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson', quiet = TRUE) |> 
  st_union() 
landcover <- rast('~/Documents/GitHub/OptimSS/data/1_raw/aligned_landcover.tif')
aedes <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/aedes_raw.geojson', quiet = TRUE) |>
  mutate(pch = 24)
anopheles <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/anopheles_raw.geojson', quiet = TRUE) |>
  mutate(pch = 21)
mosq <- rbind(aedes, anopheles)

tmap_options_reset()

(world_plot <- tm_shape(World, 
                        bbox = "FULL",
                        crs = "+proj=ortho +lat_0=10 +lon_0=120") +
    tm_polygons(fill = "white") +
    tm_shape(World %>% filter(name == 'Malaysia')) +
    tm_polygons(fill = "#ffffbf") +
    tm_style("natural", bg.color = '#e0f3f8') + 
    tm_shape(boundary %>% st_buffer(50000)) +
    tm_polygons(fill = 'red', col = 'red')  +
    tm_graticules(n.x = 20, n.y = 10, col = "black", lwd = .1, labels.show = FALSE) +
    tm_xlab("", size = 1.1) +tm_ylab("", size = 1.1)+
    #tm_layout(outer.margins = c(0.2, 0.2, 0.2, 0.2)
    NULL          )

(studyarea_plot <- ggplot() +
    geom_spatraster(data = as.factor(landcover)) +
    geom_sf(data = boundary, fill = NA, color = "black") +
    scale_fill_manual(
      values = c('#f6e8c3', '#80cdc1', '#018571', '#d8b365', 'gray'),
      labels = c("0" = "Oil palm", "1" = "Secondary forest", "2" = "Primary forest", "3" = "Other plantations", "4" = "Built-up"),
      na.value = NA,
      na.translate = FALSE,
      #guide = guide_legend(nrow = 1)
    ) +
    geom_sf(data = mosq, aes(shape = species), color = "red", fill = NA, size = 3) +
    scale_shape_manual(
      values = c("Aedes albopictus" = 24, "Anopheles balabacensis" = 21),
      labels = c("Aedes albopictus" = "Ae. albopictus", "Anopheles balabacensis" = "An. balabacensis")
    ) +
    coord_sf(xlim = c(ext(landcover)[1] - 50000, ext(landcover)[2] ),
             ylim = c(ext(landcover)[3], ext(landcover)[4])) +
    theme_minimal() +
    #labs(title = "c) Study area") +
    labs(fill = 'Land Use/Land Cover', x='', y='', shape = 'Species') +
    theme_void(base_family = 'Times New Roman', base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5), # Center the title
      #plot.background = element_rect(fill = "white", colour = 'black'),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20) # Add offset to plot boundaries
    ) +
    guides(shape = guide_legend(ncol = 1, direction = 'vertical', label.theme = element_text(face = "italic"))) +
    guides(fill = guide_legend(ncol = 1, direction = 'vertical')) +
    ggspatial::annotation_scale(location = "br", line_width = 0.5, tick_height = 0.1, text_cex = 0.6,
                                pad_y=unit(0.5,"mm")) +
    ggspatial::annotation_north_arrow(location = "tr", height = unit(2,"cm"), width = unit(2,"cm"),
                                      pad_y=unit(2,"mm"),
                                      style = ggspatial::north_arrow_nautical(line_width = 0.6, text_size = 6)) 
)

world_plot <- tmap_grob(world_plot)

(final_plot <- ggdraw() +
    draw_plot(studyarea_plot) +
    draw_plot(world_plot,
              height = .8,
              width = .8,
              x = -0.2,
              y = 0.1,
    ) +
    draw_line(
      x = c(0.2, 0.4),
      y = c(0.48, 0.8),
      color = "red", size = .5, linetype = 2
    ) +
    draw_line(
      x = c(0.2, 0.4),
      y = c(0.48, 0.2),
      color = "red", size = .5, linetype = 2
    )
)

ggsave('~/Documents/GitHub/OptimSS/data/5_figures/1_studyarea.jpg', final_plot, width = 14, height = 7, dpi = 300)
 