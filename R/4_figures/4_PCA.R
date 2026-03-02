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

# 
df_clean <- na.omit(df)
df_anopheles <- df_clean |> dplyr::filter(species == "anopheles", sample_size == 15) 
labels <- df_anopheles$scenario
# add True map as a label
labels <- c(labels, "True_Map")
metrics <- df_anopheles |> select(A_Primary_abundance, B_Secondary_abundance, C_Oil_abundance, D_Plantation_abundance, E_Built_abundance,
                                  `size for the nbinomial observations (1/overdispersion)`, `Range for field`, `Stdev for field`, `GroupRho for field`,
                                  global.global_pred_abundance, global.global_pred_contig, global.global_pred_joinent)
# add the true map metrics to the metrics data frame
metrics <- rbind(metrics, anopheles_true_df)

# Run PCA
# center = TRUE and scale. = TRUE 
# This standardizes metrics, aligning with Mahalanobis logic.
pca_res <- prcomp(metrics, center = TRUE, scale. = TRUE)

# Extract PC coordinates and calculate variance explained
pca_data <- data.frame(
  Design = labels,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3]
)

# Calculate percentage of variance explained for the axes
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

# Create a flag to separate the True Map from the Predictions for plotting
pca_data <- pca_data %>%
  mutate(Is_True_Map = ifelse(Design == "True_Map", "Yes", "No"))

# Plot
my_palette <- c(
  "True_Map" = "black",              
  "a" = '#686ea6',
  "b" = '#6da879',
  "c" = '#eba5a0',
  "d" = '#3842a1',
  "e" = '#3a9e4d',
  "f" = '#e3685f',
  "g" = '#0716a3',
  "h" = '#08a124',
  "i" = '#e61a0b'
)

create_pca_plot <- function(x_var, y_var, x_label, y_label) {
  ggplot(pca_data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                       color = Design, shape = Is_True_Map, size = Is_True_Map)) +
    geom_point(alpha = 0.8) +
    scale_shape_manual(values = c("No" = 16, "Yes" = 8)) + 
    scale_size_manual(values = c("No" = 3, "Yes" = 8)) +
    scale_color_manual(values = my_palette) + 
    labs(x = x_label, y = y_label) +
    theme_minimal() +
    guides(shape = "none", size = "none")
}

p1 <- create_pca_plot("PC1", "PC2", 
                      paste0("PC1 (", var_explained[1], "%)"), 
                      paste0("PC2 (", var_explained[2], "%)")) +
      theme(legend.position = "none")

p2 <- create_pca_plot("PC1", "PC3", 
                      paste0("PC1 (", var_explained[1], "%)"), 
                      paste0("PC3 (", var_explained[3], "%)"))

combined_plot <- p1 + p2 + 
  plot_annotation(
    title = "Multidimensional Performance of Sampling Designs",
    subtitle = "Proximity to the True Map (Black Star) indicates higher accuracy"
  )
