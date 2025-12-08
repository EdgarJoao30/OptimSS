library(tidyverse)
library(patchwork) # For combining plots
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'

# List all files in the directory
files_dir_an <- paste0(wd, '/results/0_explained_deviance/anopheles')
files_dir_ae <- paste0(wd, '/results/0_explained_deviance/aedes')
file_list_an <- list.files(files_dir_an, full.names = TRUE)
file_list_ae <- list.files(files_dir_ae, full.names = TRUE)
# Read and combine all files with sampling size extracted
an <- file_list_an %>%
  map_dfr(~ {
    read_csv(.x) 
  })
ae <- file_list_ae %>%
  map_dfr(~ {
    read_csv(.x) 
  })

ed_an <- an %>%
  filter(Effect == 'pseudo_r2_full') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_r2 = sum(Mean, na.rm = TRUE),
    sd_r2 = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  ) %>% 
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,Scenario == "b" ~ 2, Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,Scenario == "e" ~ 2,Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,Scenario == "h" ~ 2,Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,Scenario == "b" ~ 1,Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,Scenario == "e" ~ 2,Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,Scenario == "h" ~ 3,Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    )
  )

ed_ae <- ae %>%
  filter(Effect == 'pseudo_r2_full') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_r2 = sum(Mean, na.rm = TRUE),
    sd_r2 = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  ) %>% 
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,Scenario == "b" ~ 2, Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,Scenario == "e" ~ 2,Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,Scenario == "h" ~ 2,Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,Scenario == "b" ~ 1,Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,Scenario == "e" ~ 2,Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,Scenario == "h" ~ 3,Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    )
  )


an_summ <- ed_an %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(mean_r2 = mean(mean_r2), 
            gradient1 = max(gradient1), gradient2 = max(gradient2)) %>% 
  group_by(sample_size) %>% 
  slice_max(order_by = mean_r2, n = 1, with_ties = FALSE) %>% 
  mutate(y = 0.32)

ae_summ <- ed_ae %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(mean_r2 = mean(mean_r2), 
            gradient1 = max(gradient1), gradient2 = max(gradient2)) %>% 
  group_by(sample_size) %>% 
  slice_max(order_by = mean_r2, n = 1, with_ties = FALSE) %>% 
  mutate(y = 0.58)


(g_an <- ggplot(data = ed_an, 
                aes(x = as.factor(sample_size), y = mean_r2)) +
    # geom_rect(aes(colour = back_fill),
    #           xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = NA, linetype = 2, size = .5)+
    scale_colour_identity()+
    geom_boxplot() +
    facet_grid(gradient2~gradient1, scales = "fixed", 
               labeller = labeller(
                 gradient1 = c("1" = "Uniform", "2" = "Stratified", "3" = "Random"),
                 gradient2 = c("1" = "Fixed", "2" = "Mixed", "3" = "Variable")
               )) +
    
    labs(
      title = expression(italic("(a) Anopheles balabacensis")),
      x = "Sample Size",
      y = "pseudo R2"
    ) +
    theme_minimal(base_family = 'Times New Roman', base_size = 12) +
    # theme(legend.position = 'bottom',
    #       #strip.background = element_rect(fill = "white", color = "black"),
    #       strip.text = element_text(color = "black")) +
    # scale_y_log10(breaks = c(1e0, 1e-1, 1e-2, 1e-3)) + 
    geom_text(
      data = an_summ,
      aes(x = factor(sample_size), y = y ),
      label = "*",
      size = 8,
      color = "red"
    ) +
    ylim(0, .35)+
    NULL
)

(g_ae <- ggplot(data = ed_ae, 
                aes(x = as.factor(sample_size), y = mean_r2)) +
    # geom_rect(aes(colour = back_fill),
    #           xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = NA, linetype = 2, size = .5)+
    scale_colour_identity()+
    geom_boxplot() +
    facet_grid(gradient2~gradient1, scales = "fixed", 
               labeller = labeller(
                 gradient1 = c("1" = "Uniform", "2" = "Stratified", "3" = "Random"),
                 gradient2 = c("1" = "Fixed", "2" = "Mixed", "3" = "Variable")
               )) +
    
    labs(
      title = expression(italic("(b) Aedes albopictus")),
      x = "Sample Size",
      y = "pseudo R2"
    ) +
    theme_minimal(base_family = 'Times New Roman', base_size = 12) +
    # theme(legend.position = 'bottom',
    #       #strip.background = element_rect(fill = "white", color = "black"),
    #       strip.text = element_text(color = "black")) +
    # scale_y_log10(breaks = c(1e0, 1e-1, 1e-2, 1e-3)) + 
    geom_text(
      data = ae_summ,
      aes(x = factor(sample_size), y = y ),
      label = "*",
      size = 8,
      color = "red"
    ) +
    ylim(0, .65)+
    NULL
)

(final_plot <- g_an | g_ae + plot_layout(axis_titles = "collect")) 

ggsave(paste0(wd, '/figures/r2.jpg'), final_plot, dpi = 300, width = 7, height = 7)
