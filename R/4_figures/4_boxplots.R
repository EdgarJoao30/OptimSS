packs <- c("tidyverse", "RColorBrewer", "patchwork")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

# load files
anopheles_folder <- "~/Documents/GitHub/OptimSS/data/4_results/results/anopheles/"
aedes_folder <- "~/Documents/GitHub/OptimSS/data/4_results/results/aedes/"
anopheles_files <- list.files(anopheles_folder, pattern = "anopheles", full.names = TRUE)
aedes_files <- list.files(aedes_folder, pattern = "aedes", full.names = TRUE)
anopheles_data <- lapply(anopheles_files, read_csv) |> bind_rows()
aedes_data <- lapply(aedes_files, read_csv) |> bind_rows()
df <- rbind(anopheles_data, aedes_data) 
anopheles_true <- read.csv("~/Documents/GitHub/OptimSS/data/1_raw/anopheles_parameters.csv") 
#
anopheles_true_df <- anopheles_true |>
     dplyr::select(ID, mean_abundance) |>
     spread(key = ID, value = mean_abundance) 
colnames(anopheles_true_df) <- c("A_Primary_abundance", "B_Secondary_abundance", "C_Oil_abundance", "D_Plantation_abundance", "E_Built_abundance")
anopheles_true_df$`size for the nbinomial observations (1/overdispersion)` <- 10
anopheles_true_df$`Range for field` <- 3000
anopheles_true_df$`Stdev for field` <- 0.5
anopheles_true_df$`GroupRho for field` <- 0.9
anopheles_true_df$global.global_pred_abundance <- unique(anopheles_data$global.global_sim_abundance)
anopheles_true_df$global.global_pred_contig <- unique(anopheles_data$global.global_sim_contig)
anopheles_true_df$global.global_pred_joinent <- unique(anopheles_data$global.global_sim_joinent)

colnames(df)

(g_an <- ggplot(data = df |> filter(species == "anopheles"), 
                aes(x = as.factor(scenario), y = global.global_frac_ks)) +
    scale_colour_identity() +
    geom_boxplot() +
    facet_grid(sample_size ~ ., scales = "fixed") +
    labs(
      title = expression(italic("(a) Anopheles balabacensis")),
      x = "Scenario",
      y = "RMSE of global abundance"
    ) +
    theme_minimal(base_family = 'Times New Roman', base_size = 12) +
    ylim(0, .5) +
    NULL
)
