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

colnames(anopheles_true_df) <- c(
  "A_Primary_abundance",
  "B_Secondary_abundance",
  "C_Oil_abundance",
  "D_Plantation_abundance",
  "E_Built_abundance"
)

anopheles_true_df <- anopheles_true_df |>
  mutate(
    `size for the nbinomial observations (1/overdispersion)` = 10,
    `Range for field` = 3000,
    `Stdev for field` = 0.5,
    `GroupRho for field` = 0.9,
    global.global_abundance_rmse = 0,
    global.global_contig_rmse = 0,
    global.global_joinent_rmse = 0
  )

# 
df_clean <- na.omit(df)
# remove rows with outliers
df_clean <- df_clean |> filter(focal.focal_abundance_rmse < 15)
# remove rows with outliers in Range for field
df_clean <- df_clean |> filter(`Range for field` < 5000)
df_anopheles <- df_clean |> dplyr::filter(species == "anopheles", sample_size == 15) 
labels <- df_anopheles$scenario
# add True map as a label
labels <- c(labels, "True_Map")
# Global metrics
metrics_global <- df_anopheles |> dplyr::select(A_Primary_abundance, B_Secondary_abundance, C_Oil_abundance, D_Plantation_abundance, E_Built_abundance,
                                  `size for the nbinomial observations (1/overdispersion)`, `Range for field`, `Stdev for field`, `GroupRho for field`,
                                  global.global_abundance_rmse, global.global_contig_rmse, global.global_joinent_rmse)

metrics_global$A_Primary_abundance <- abs(metrics_global$A_Primary_abundance - anopheles_true_df$A_Primary_abundance)
metrics_global$B_Secondary_abundance <- abs(metrics_global$B_Secondary_abundance - anopheles_true_df$B_Secondary_abundance)
metrics_global$C_Oil_abundance <- abs(metrics_global$C_Oil_abundance - anopheles_true_df$C_Oil_abundance)
metrics_global$D_Plantation_abundance <- abs(metrics_global$D_Plantation_abundance - anopheles_true_df$D_Plantation_abundance)
metrics_global$E_Built_abundance <- abs(metrics_global$E_Built_abundance - anopheles_true_df$E_Built_abundance)
metrics_global$`size for the nbinomial observations (1/overdispersion)` <- abs(metrics_global$`size for the nbinomial observations (1/overdispersion)` - anopheles_true_df$`size for the nbinomial observations (1/overdispersion)`)
metrics_global$`Range for field` <- abs(metrics_global$`Range for field` - anopheles_true_df$`Range for field`)
metrics_global$`Stdev for field` <- abs(metrics_global$`Stdev for field` - anopheles_true_df$`Stdev for field`)
# Local and focal metrics
metrics_local <- df_anopheles |> dplyr::select(local.local_abundance_rmse, 
                                               focal.focal_abundance_rmse,
                                               focal.focal_contig_rmse,
                                               focal.focal_joinent_rmse, 
                                               focal.focal_classfreq_rmse)

# true map row for global metrics
true_map_row_global <- data.frame(
  A_Primary_abundance = 0,
  B_Secondary_abundance = 0,
  C_Oil_abundance = 0,
  D_Plantation_abundance = 0,
  E_Built_abundance = 0,
  "size for the nbinomial observations (1/overdispersion)" = 0,
  "Range for field" = 0,
  "Stdev for field" = 0,
  "GroupRho for field" = 0,
  global.global_abundance_rmse = 0,
  global.global_contig_rmse = 0,
  global.global_joinent_rmse = 0
)
colnames(true_map_row_global) <- colnames(metrics_global)

# true map row for local and focal metrics
true_map_row_local <- data.frame(
  local.local_abundance_rmse = 0,
  focal.focal_abundance_rmse = 0,
  focal.focal_contig_rmse = 0,
  focal.focal_joinent_rmse = 0,
  focal.focal_classfreq_rmse = 0
)
colnames(true_map_row_local) <- colnames(metrics_local)

metrics_global <- rbind(metrics_global, true_map_row_global)
metrics_local <- rbind(metrics_local, true_map_row_local)

# add the true map metrics to the metrics data frame
# metrics <- rbind(metrics, anopheles_true_df)
# check correlation between metrics
cor_matrix <- cor(metrics_local)
print(cor_matrix)
# plot correlation matrix
cor_matrix_long <- as.data.frame(as.table(cor_matrix))
ggplot(cor_matrix_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()
# drop columns that are highly correlated with each other (correlation > 0.9)
# metrics <- metrics |> dplyr::select(-global.global_pred_abundance, -global.global_pred_contig, -global.global_pred_joinent)
# metrics <- metrics |> dplyr::select(`size for the nbinomial observations (1/overdispersion)`:global.global_pred_joinent)
# Run PCA
# center = TRUE and scale. = TRUE 
# This standardizes metrics, aligning with Mahalanobis logic.
pca_res_global <- prcomp(metrics_global, center = TRUE, scale. = TRUE)
pca_res_local <- prcomp(metrics_local, center = TRUE, scale. = TRUE)
# Extract PC coordinates and calculate variance explained
pca_data_global <- data.frame(
  Design = labels,
  PC1 = pca_res_global$x[, 1],
  PC2 = pca_res_global$x[, 2],
  PC3 = pca_res_global$x[, 3]
)
pca_data_local <- data.frame(
  Design = labels,
  PC1 = pca_res_local$x[, 1],
  PC2 = pca_res_local$x[, 2],
  PC3 = pca_res_local$x[, 3]
)

# Calculate percentage of variance explained for the axes
var_explained_global <- round(100 * pca_res_global$sdev^2 / sum(pca_res_global$sdev^2), 1)
var_explained_local <- round(100 * pca_res_local$sdev^2 / sum(pca_res_local$sdev^2), 1)

# Create a flag to separate the True Map from the Predictions for plotting
pca_data_global <- pca_data_global %>%
  mutate(Is_True_Map = ifelse(Design == "True_Map", "Yes", "No"))
pca_data_local <- pca_data_local %>%
  mutate(Is_True_Map = ifelse(Design == "True_Map", "Yes", "No"))
# Plot
my_palette <- c(
  "True_Map" = "black",              
  "a" = '#686ea6',
  "b" = '#3842a1',
  "c" = '#0716a3',
  "d" = '#eba5a0',
  "e" = '#e3685f',
  "f" = '#e61a0b',
  "g" = '#6da879',
  "h" = '#3a9e4d',
  "i" = '#08a124'
)
nrow(pca_data_global)
pca_data <- cbind(pca_data_global |> dplyr::select(Design, Is_True_Map, PC1) |> rename(PC1_global = PC1), 
                  pca_data_local |> dplyr::select(PC1_local = PC1))

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

p1 <- create_pca_plot("PC1_global", "PC1_local", 
                      paste0("PC1 Global (", var_explained_global[1], "%)"), 
                      paste0("PC1 Local (", var_explained_local[1], "%)")) +
      # theme(legend.position = "none") +
      NULL

p2 <- create_pca_plot("PC3", "PC2", 
                      paste0("PC3 (", var_explained[3], "%)"),
                      paste0("PC2 (", var_explained[2], "%)"))

combined_plot <- p1 + p2 + 
  plot_annotation(
    title = "Multidimensional Performance of Sampling Designs",
    subtitle = "Proximity to the True Map (Black Star) indicates higher accuracy"
  )
