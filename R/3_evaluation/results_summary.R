library(tidyverse)
library(patchwork) # For combining plots
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'

# List all files in the directory
files_dir_ae <- paste0(wd, '/results/aedes')
files_dir_an <- paste0(wd, '/results/anopheles')

file_list_ae <- list.files(files_dir_ae, full.names = TRUE)
file_list_an <- list.files(files_dir_an, full.names = TRUE)

# Read and combine all files with sampling size extracted
ae <- file_list_ae %>%
  map_dfr(~ {
    read_csv(.x) 
  })

an <- file_list_an %>%
  map_dfr(~ {
    read_csv(.x) 
  })

ae <- ae %>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    )
  )

an <- an %>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    )
  )

# Aggregate the coefficients' mean and sd across iterations
# 0.7351552, 2.588883, 1.775618, 1.612111, 1.156367
summary_coeffs <- an %>%
  group_by(sample_size, Scenario, Effect) %>%
  summarise(
    mean = mean(Mean, na.rm = TRUE),
    sd = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )

# Create an index for each Scenario using a loss function


coeff_loss_an <- an %>%
  filter(Effect %in% c("Intercept", 
                       "B_Secondary", 
                       "C_Oil", "D_Plantation", 
                       "E_Built", "size for the nbinomial observations (1/overdispersion)",
                       "Range for field",
                       "GroupRho for field")) %>%
  mutate(
    true_value = case_when(
      Effect == "Intercept" ~ 0.7351552,
      Effect == "B_Secondary" ~ 2.588883,
      Effect == "C_Oil" ~ 1.775618,
      Effect == "D_Plantation" ~ 1.612111,
      Effect == "E_Built" ~ 1.156367,
      Effect == "size for the nbinomial observations (1/overdispersion)" ~ 10,
      Effect == "Range for field" ~ 3000,
      Effect == "GroupRho for field" ~ 0.935,
      TRUE ~ NA_real_
    ),
    coefficient_loss = (Mean - true_value)^2
  ) %>%
  group_by(Effect) %>%
  mutate(
    coefficient_loss_norm = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE))) %>% 
  group_by(sample_size, Scenario, iteration) %>% 
  summarise(
    coefficient_loss = sum(coefficient_loss_norm, na.rm = TRUE),
    .groups = "drop"
  ) 
  # group_by(sample_size, Scenario) %>%
  # #mutate(coefficient_loss_norm = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE))) %>%
  # summarise(
  #   coefficient_loss = mean(coefficient_loss, na.rm = TRUE),
  #   .groups = "drop"
  # )


rmse_an <- an %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_rmse = sum(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )
index_an <- coeff_loss_an
index_an$rmse <- rmse_an$mean_rmse

index_an_ns <- index_an %>%
  mutate(
    norm_rmse = (rmse - min(rmse, na.rm = TRUE)) / (max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE)),
    norm_coefficient_loss = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE)),
    index = (0.5 + (norm_rmse + norm_coefficient_loss)) * sqrt(sample_size)
  )

index_an_ns_s <- index_an_ns %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(index = mean(index))

index_an_ns_s_2 <- index_an_ns_s %>%  group_by(sample_size) %>% 
  slice_min(order_by = index, n = 1, with_ties = FALSE)

index_an <- index_an %>%
  group_by(sample_size, Scenario) %>% 
  summarise(rmse = mean(rmse),
            coefficient_loss = mean(coefficient_loss)) %>% 
  mutate(
    norm_rmse = (rmse - min(rmse, na.rm = TRUE)) / (max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE)),
    norm_coefficient_loss = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE)),
    index = (0.5 + (norm_rmse + norm_coefficient_loss)) * sqrt(sample_size)
  )

index_an_2 <- index_an %>%  group_by(sample_size) %>% 
  slice_min(order_by = index, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    y = 40
  )
summary_rmse <- an %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_rmse = mean(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )

rmse_an <- summary_rmse %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(rmse = mean(mean_rmse, na.rm = TRUE))

rmse_an_2 <- rmse_an %>%  group_by(sample_size) %>% 
  slice_min(order_by = rmse, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    y = 40
  )

# Aedes
summary_rmse_ae <- ae %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_rmse = mean(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )

coeff_loss_ae <- ae %>%
  filter(Effect %in% c("Intercept", 
                       "B_Secondary", 
                       "C_Oil", "D_Plantation", 
                       "E_Built", "size for the nbinomial observations (1/overdispersion)",
                       "Range for field",
                       "GroupRho for field")) %>%
  mutate(
    true_value = case_when(
      Effect == "Intercept" ~ -0.4310554,
      Effect == "B_Secondary" ~ 0.50800013,
      Effect == "C_Oil" ~ 0.06971481,
      Effect == "D_Plantation" ~ 0.30783431,
      Effect == "E_Built" ~ 3.30660488,
      Effect == "size for the nbinomial observations (1/overdispersion)" ~ 10,
      Effect == "Range for field" ~ 1000,
      Effect == "GroupRho for field" ~ 0.9,
      TRUE ~ NA_real_
    ),
    coefficient_loss = (Mean - true_value)^2
  ) %>%
  group_by(Effect) %>%
  mutate(
    coefficient_loss_norm = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE))) %>% 
  group_by(sample_size, Scenario) %>%
  summarise(
    coefficient_loss = sum(coefficient_loss_norm, na.rm = TRUE),
    .groups = "drop"
  ) 

rmse_ae <- ae %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario) %>%
  summarise(
    mean_rmse = mean(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )
index_ae <- coeff_loss_ae
index_ae$rmse <- rmse_ae$mean_rmse

index_ae <- index_ae %>%
  #group_by(sample_size, Scenario) %>% 
  mutate(
    norm_rmse = (rmse - min(rmse, na.rm = TRUE)) / (max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE)),
    norm_coefficient_loss = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE)),
    index = (1 + (norm_rmse + norm_coefficient_loss)) * sqrt(sample_size )
  )






rmse_ae <- summary_rmse_ae %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(rmse = mean(mean_rmse, na.rm = TRUE))

rmse_ae_2 <- rmse_ae %>%  group_by(sample_size) %>% 
  slice_min(order_by = rmse, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    y = 10
  )

index_ae_2 <- index_ae %>%  group_by(sample_size) %>% 
  slice_min(order_by = index, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(
    gradient1 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 2,
      Scenario == "c" ~ 3,
      Scenario == "d" ~ 1,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 3,
      Scenario == "g" ~ 1,
      Scenario == "h" ~ 2,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    gradient2 = case_when(
      Scenario == "a" ~ 1,
      Scenario == "b" ~ 1,
      Scenario == "c" ~ 1,
      Scenario == "d" ~ 2,
      Scenario == "e" ~ 2,
      Scenario == "f" ~ 2,
      Scenario == "g" ~ 3,
      Scenario == "h" ~ 3,
      Scenario == "i" ~ 3,
      TRUE ~ NA_real_
    ),
    y = 10
  )
mean_rmse <- summary_rmse %>% group_by(Scenario, sample_size) %>% summarise(mean_rmse = median(mean_rmse))

# Calculate the global min and max for the fill scale
global_min <- min(mean_rmse$mean_rmse, na.rm = TRUE)
global_max <- max(mean_rmse$mean_rmse, na.rm = TRUE)

# Define the true values for all coefficients
params_an <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367, 3000, 0.935, 10)
params_ae <- c(-0.4310554, 0.50800013, 0.06971481, 0.30783431, 3.30660488, 1000, 0.9, 10)

true_values_an <- tibble(
  Effect = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "Range for field", "GroupRho for field", "size for the nbinomial observations (1/overdispersion)"),
  true_value = params_an
) %>% 
  mutate(Effect = factor(Effect, levels = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "GroupRho for field", "Range for field", "size for the nbinomial observations (1/overdispersion)")),
         order_effect = case_when(
           Effect == "Intercept" ~ 1,
           Effect == "B_Secondary" ~ 2,
           Effect == "C_Oil" ~ 3,
           Effect == "D_Plantation" ~ 4,
           Effect == "E_Built" ~ 5,
           Effect == "GroupRho for field" ~ 6,
           Effect == "Range for field" ~ 7, 
           Effect == "size for the nbinomial observations (1/overdispersion)" ~ 8,
           TRUE ~ NA_real_
         ))

true_values_ae <- tibble(
  Effect = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "Range for field", "GroupRho for field", "size for the nbinomial observations (1/overdispersion)"),
  true_value = params_ae
) %>% 
  mutate(Effect = factor(Effect, levels = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "GroupRho for field", "Range for field", "size for the nbinomial observations (1/overdispersion)")),
         order_effect = case_when(
           Effect == "Intercept" ~ 1,
           Effect == "B_Secondary" ~ 2,
           Effect == "C_Oil" ~ 3,
           Effect == "D_Plantation" ~ 4,
           Effect == "E_Built" ~ 5,
           Effect == "GroupRho for field" ~ 6,
           Effect == "Range for field" ~ 7, 
           Effect == "size for the nbinomial observations (1/overdispersion)" ~ 8,
           TRUE ~ NA_real_
         ))

# Merge the true values with the summary data
summary_coeffs <- summary_coeffs %>%
  inner_join(true_values, by = "Effect")

# RMSE plot
(g_an <- ggplot(an %>% filter(Effect == 'rmse'), 
  aes(x = factor(sample_size), y = Mean)) +
  geom_boxplot() +
  ylim(NA, 40) + # 50 for anopheles, 10 for aedes
  facet_grid(gradient2~gradient1, scales = "fixed", 
    labeller = labeller(
    gradient1 = c("1" = "Uniform", "2" = "Stratified", "3" = "Random"),
    gradient2 = c("1" = "Fixed", "2" = "Mixed", "3" = "Variable")
    )) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(summary_rmse$mean_rmse)) +
  labs(
    title = expression("A. " ~ italic("Anopheles balabacensis")),
    x = "Sample Size",
    y = "Overall RMSE",
    fill = "Overall RMSE"
  ) +
  theme_minimal() +
  labs(
    #caption = "Columns: Spatial structure \nRows: Temporal flexibility"
  ) +
  theme(legend.position = 'bottom') +
  geom_text(
    data = index_an_2,
    aes(x = factor(sample_size), y = y -2), # adjust y as needed
    label = "*",
    size = 8,
    color = "red"
  ) +
    geom_text(
      data = index_an_2 %>% slice_min(order_by = index, n = 1, with_ties = FALSE),
      aes(x = factor(sample_size), y = y -1), # adjust y as needed
      label = "O",
      size = 8,
      color = "red"
    )
)

(g_ae <- ggplot(ae %>% filter(Effect == 'rmse'), 
                aes(x = factor(sample_size), y = Mean,)) +
    geom_boxplot() +
    ylim(NA, 10) + # 50 for anopheles, 10 for aedes
    facet_grid(gradient2~gradient1, scales = "fixed", 
               labeller = labeller(
                 gradient1 = c("1" = "Uniform", "2" = "Stratified", "3" = "Random"),
                 gradient2 = c("1" = "Fixed", "2" = "Mixed", "3" = "Variable")
               )) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(summary_rmse_ae$mean_rmse)) +
    labs(
      title = expression("B. " ~ italic("Aedes albopictus")),
      x = "Sample Size",
      y = "Overall RMSE",
      fill = "Overall RMSE"
    ) +
    theme_minimal() +
    labs(
      #caption = "Columns: Spatial structure \nRows: Temporal flexibility"
    ) +
    theme(legend.position = 'bottom')+
    geom_text(
      data = index_ae_2,
      aes(x = factor(sample_size), y = y - .58), # adjust y as needed
      label = "*",
      size = 8,
      color = "red"
    )+
    geom_text(
      data = index_ae_2 %>% slice_min(order_by = index, n = 1, with_ties = FALSE),
      aes(x = factor(sample_size), y = y -.3), # adjust y as needed
      label = "O",
      size = 8,
      color = "red"
    )
)

(rmse_plot <- g_an | g_ae)

# ggsave(paste0(wd, '/figures/rmse3.jpg'), rmse_plot, dpi = 300, width = 7, height = 7)

plots <- combined_data %>%
  filter(!(Effect %in% c("Stdev for field", "rmse"))) %>% 
  mutate(Effect = factor(Effect, levels = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "GroupRho for field", "Range for field", "size for the nbinomial observations (1/overdispersion)")),
         order_effect = case_when(
           Effect == "Intercept" ~ 1,
           Effect == "B_Secondary" ~ 2,
           Effect == "C_Oil" ~ 3,
           Effect == "D_Plantation" ~ 4,
           Effect == "E_Built" ~ 5,
           Effect == "GroupRho for field" ~ 6,
           Effect == "Range for field" ~ 7,
           Effect == "size for the nbinomial observations (1/overdispersion)" ~ 8,
           TRUE ~ NA_real_
         )) %>% 
  split(.$Scenario) %>%
  map(~ {
    scenario_rmse <- mean_rmse %>% filter(Scenario == unique(.x$Scenario))
    ggplot(.x,  
           aes(x = as.factor(sample_size), y = Mean, fill = scenario_rmse$mean_rmse[match(.x$sample_size, scenario_rmse$sample_size)])) +
      geom_boxplot() +
      facet_wrap(~order_effect, scale = "fixed", nrow = 1, labeller = labeller(order_effect = c(
        "1" = "Primary",
        "2" = "Secondary",
        "3" = "Oil palm",
        "4" = "Other plantation",
        "5" = "Built-up",
        "6" = "Rho",
        "7" = "Spatial range",
        "8" = "Overdispersion"
      ))) +
      geom_hline(data = true_values2, aes(yintercept = true_value), color = "red") +
      labs(
        title = paste("Scenario:", unique(.x$Scenario)),
        x = "",
        y = "",
        fill = "Overall RMSE"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
        strip.text = element_text(angle = 65, hjust = 1,),   # Rotate facet titles for better readability
        legend.position = "none") +
      scale_y_log10() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = median(mean_rmse$mean_rmse, na.rm = TRUE), 
                           limits = c(global_min, global_max))
  })

# Combine all plots into a grid
combined_plot <- wrap_plots(plots, ncol = 3)
combined_plot

# ggsave(paste0(wd, '/figures/coeffs_rmse.jpg'), combined_plot, dpi = 300, width = 16, height = 9)



# Plots per effect

effect_plots <- ae %>%
  filter(!(Effect %in% c("Stdev for field", "rmse"))) %>% 
  mutate(Effect = factor(Effect, levels = c("Intercept", "B_Secondary", "C_Oil", "D_Plantation", "E_Built", "GroupRho for field", "Range for field", "size for the nbinomial observations (1/overdispersion)")),
         order_effect = case_when(
           Effect == "Intercept" ~ 1,
           Effect == "B_Secondary" ~ 2,
           Effect == "C_Oil" ~ 3,
           Effect == "D_Plantation" ~ 4,
           Effect == "E_Built" ~ 5,
           Effect == "GroupRho for field" ~ 6,
           Effect == "Range for field" ~ 7,
           Effect == "size for the nbinomial observations (1/overdispersion)" ~ 8,
           TRUE ~ NA_real_
         )) %>% 
      mutate(Effect = recode(Effect,
           "Intercept" = "Primary forest",
           "B_Secondary" = "Secondary forest",
           "C_Oil" = "Oil palm",
           "D_Plantation" = "Other plantations",
           "E_Built" = "Built-up",
           "GroupRho for field" = "Temporal correlation",
           "Range for field" = "Spatial range",
           "size for the nbinomial observations (1/overdispersion)" = "Overdispersion")) %>% 
  split(.$order_effect) %>%
  map(~ {
    coeff_value <- true_values_ae %>% filter(order_effect == unique(.x$order_effect))
    p <- ggplot(.x,  
           aes(x = as.factor(sample_size), y = Mean)) +
      geom_violin() +
      facet_grid(gradient2~gradient1, scale = "fixed", 
                 labeller = labeller(
                   gradient1 = c("1" = "Uniform", "2" = "Stratified", "3" = "Random"),
                   gradient2 = c("1" = "Fixed", "2" = "Mixed", "3" = "Variable")
                 )) +
      geom_hline(aes(yintercept = coeff_value$true_value), color = "red") +
      labs(
        title = unique(.x$Effect),
        x = "",
        y = ""
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 12),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        legend.position = "none")
    if (unique(.x$Effect) != "Primary forest") {
      p <- p + scale_y_log10()
    }
    p
  })

# Combine all plots into a grid
effect_plots <- wrap_plots(effect_plots,nrow = 4)
effect_plots

# ggsave(paste0(wd, '/figures/coeffs_ae.jpg'), effect_plots, dpi = 300, width = 7, height = 10)
