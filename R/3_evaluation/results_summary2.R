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

# Create an index for each Scenario using a loss function

### Anopheles 
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


rmse_an <- an %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_rmse = sum(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )

loss_an <- coeff_loss_an
loss_an$rmse <- rmse_an$mean_rmse

loss_an2 <- loss_an %>%
  mutate(
    norm_rmse = (rmse - min(rmse, na.rm = TRUE)) / (max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE)),
    norm_coefficient_loss = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE)),
    index = ( (norm_rmse + norm_coefficient_loss)) #* sqrt(sample_size)
  ) %>% 
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
  ) %>% 
  mutate(back_fill = ifelse(Scenario == 'b', 'red', 'white'))

loss_an_summ <- loss_an2 %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(loss_mean = mean(index), 
            loss_sum = sum(index),
            loss_0.025 = quantile(index, 0.025),
            loss_0.975 = quantile(index, 0.975),
            rmse_mean = mean(rmse), 
            rmse_0.025 = quantile(rmse, 0.025),
            rmse_0.975 = quantile(rmse, 0.975),
            pl_mean = mean(coefficient_loss), 
            pl_0.025 = quantile(coefficient_loss, 0.025),
            pl_0.975 = quantile(coefficient_loss, 0.975),
            gradient1 = max(gradient1), gradient2 = max(gradient2))

# write_csv(loss_an_summ, paste0(wd, '/results/loss_an_summ.csv'))

loss_an_summ2 <- loss_an_summ %>%  group_by(sample_size) %>% 
  slice_min(order_by = loss_mean, n = 1, with_ties = FALSE) %>% 
  mutate(y = 0.8)

(g_an <- ggplot(data = loss_an2, 
                aes(x = as.factor(sample_size), y = index)) +
    geom_rect(aes(colour = back_fill),
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = NA, linetype = 2, size = .5)+
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
      y = "Overall Loss"
    ) +
    theme_minimal(base_family = 'Times New Roman', base_size = 12) +
    theme(legend.position = 'bottom',
          #strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(color = "black")) +
    scale_y_log10(breaks = c(1e0, 1e-1, 1e-2, 1e-3)) + 
    geom_text(
      data = loss_an_summ2,
      aes(x = factor(sample_size), y = y ),
      label = "*",
      size = 8,
      color = "red"
    ) +
    NULL
)

### Aedes

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
  group_by(sample_size, Scenario, iteration) %>% 
  summarise(
    coefficient_loss = sum(coefficient_loss_norm, na.rm = TRUE),
    .groups = "drop"
  ) 

rmse_ae <- ae %>%
  filter(Effect == 'rmse') %>% 
  group_by(sample_size, Scenario, iteration) %>%
  summarise(
    mean_rmse = mean(Mean, na.rm = TRUE),
    sd_rmse = sqrt(sum(SD^2 + (Mean - mean(Mean, na.rm = TRUE))^2, na.rm = TRUE) / n()),
    .groups = "drop"
  )

loss_ae <- coeff_loss_ae
loss_ae$rmse <- rmse_ae$mean_rmse

loss_ae2 <- loss_ae %>%
  mutate(
    norm_rmse = (rmse - min(rmse, na.rm = TRUE)) / (max(rmse, na.rm = TRUE) - min(rmse, na.rm = TRUE)),
    norm_coefficient_loss = (coefficient_loss - min(coefficient_loss, na.rm = TRUE)) / (max(coefficient_loss, na.rm = TRUE) - min(coefficient_loss, na.rm = TRUE)),
    index = ( (norm_rmse + norm_coefficient_loss)) #* sqrt(sample_size)
  ) %>% 
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
  ) %>% 
  mutate(back_fill = ifelse(Scenario == 'b', 'red', 'white'))

loss_ae_summ <- loss_ae2 %>% 
  group_by(sample_size, Scenario) %>% 
  summarise(loss_mean = mean(index), 
            loss_sum = sum(index),
            loss_0.025 = quantile(index, 0.025),
            loss_0.975 = quantile(index, 0.975),
            rmse_mean = mean(rmse), 
            rmse_0.025 = quantile(rmse, 0.025),
            rmse_0.975 = quantile(rmse, 0.975),
            pl_mean = mean(coefficient_loss), 
            pl_0.025 = quantile(coefficient_loss, 0.025),
            pl_0.975 = quantile(coefficient_loss, 0.975),
            gradient1 = max(gradient1), gradient2 = max(gradient2))


summ1 <- loss_an_summ %>% group_by(Scenario) %>% summarize(loss_sum = sum(loss_sum))
summ2 <- loss_ae_summ %>% group_by(Scenario) %>% summarize(loss_sum = sum(loss_sum))

summ1$loss_sum_ae <- summ2$loss_sum
summ1$sum <- summ1$loss_sum + summ1$loss_sum_ae

summ1 <- summ1 %>% 
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
    )) %>% 
  mutate(sum_fill = ifelse(Scenario == 'b', 'red', 'white'))
# write_csv(loss_ae_summ, paste0(wd, '/results/loss_ae_summ.csv'))

loss_ae_summ2 <- loss_ae_summ %>%  group_by(sample_size) %>% 
  slice_min(order_by = loss_mean, n = 1, with_ties = FALSE) %>% 
  mutate(y = 0.8)

(g_ae <- ggplot(loss_ae2, 
                aes(x = factor(sample_size), y = index)) +
    geom_rect(aes(colour = back_fill),
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = NA, linetype = 2, size = .5)+
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
      y = " Overall Loss"
    ) +
    theme_minimal(base_family = 'Times New Roman', base_size = 12) +
    theme(legend.position = 'bottom') +
    scale_y_log10(breaks = c(1e0, 1e-1, 1e-2, 1e-3)) + 
    geom_text(
      data = loss_ae_summ2,
      aes(x = factor(sample_size), y = y ),
      label = "*",
      size = 8,
      color = "red"
    ) 
)

(final_plot <- g_an | g_ae + plot_layout(axis_titles = "collect")) 

# ggsave(paste0(wd, '/figures/loss2.jpg'), final_plot, dpi = 300, width = 7, height = 7)
