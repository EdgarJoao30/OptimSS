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

# Keep failed samples (NA estimates) and handle them explicitly downstream.
df_clean <- df
# remove rows with outliers
df_clean <- df_clean |> filter(is.na(focal.focal_abundance_rmse) | focal.focal_abundance_rmse < 15)
# remove rows with outliers in Range for field
df_clean <- df_clean |> filter(is.na(`Range for field`) | `Range for field` < 5000)
df_anopheles <- df_clean |> dplyr::filter(species == "anopheles", sample_size == 15) 
labels <- df_anopheles$scenario
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
metrics_global$`GroupRho for field` <- abs(metrics_global$`GroupRho for field` - anopheles_true_df$`GroupRho for field`)
# Local and focal metrics
metrics_local <- df_anopheles |> dplyr::select(local.local_abundance_rmse, 
                                               focal.focal_abundance_rmse,
                                               focal.focal_contig_rmse,
                                               focal.focal_joinent_rmse, 
                                               focal.focal_classfreq_rmse)

# Missing estimates represent sampling limitations; track and penalize them.
metrics_global_na_count <- rowSums(is.na(metrics_global))

impute_with_penalty <- function(df_metrics, multiplier = 1.10) {
  df_out <- df_metrics
  for (j in seq_along(df_out)) {
    col_vals <- df_out[[j]]
    col_max <- suppressWarnings(max(col_vals, na.rm = TRUE))
    if (!is.finite(col_max)) {
      col_max <- 1
    }
    penalty_val <- col_max * multiplier
    col_vals[is.na(col_vals)] <- penalty_val
    df_out[[j]] <- col_vals
  }
  df_out
}

metrics_global <- impute_with_penalty(metrics_global)

metrics_global$failed_metric_count <- metrics_global_na_count

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
  global.global_joinent_rmse = 0,
  failed_metric_count = 0
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

# Run PCA
# center = TRUE and scale. = TRUE 
# This standardizes metrics, aligning with Mahalanobis logic.
pca_res_global <- prcomp(metrics_global, center = TRUE, scale. = TRUE)
pca_res_local <- prcomp(metrics_local, center = TRUE, scale. = TRUE)

# Project the true map as a supplementary point so it does not affect PCA rotation.
true_global_scaled <- scale(true_map_row_global, center = pca_res_global$center, scale = pca_res_global$scale)
true_local_scaled <- scale(true_map_row_local, center = pca_res_local$center, scale = pca_res_local$scale)
true_global_scores <- true_global_scaled %*% pca_res_global$rotation
true_local_scores <- true_local_scaled %*% pca_res_local$rotation

# Extract PC coordinates and calculate variance explained
pca_data_global <- data.frame(
  Design = labels,
  PC1 = pca_res_global$x[, 1],
  PC2 = pca_res_global$x[, 2],
  PC3 = pca_res_global$x[, 3]
)
pca_data_global <- rbind(
  pca_data_global,
  data.frame(
    Design = "True_Map",
    PC1 = true_global_scores[1, 1],
    PC2 = true_global_scores[1, 2],
    PC3 = true_global_scores[1, 3]
  )
)

pca_data_local <- data.frame(
  Design = labels,
  PC1 = pca_res_local$x[, 1],
  PC2 = pca_res_local$x[, 2],
  PC3 = pca_res_local$x[, 3]
)
pca_data_local <- rbind(
  pca_data_local,
  data.frame(
    Design = "True_Map",
    PC1 = true_local_scores[1, 1],
    PC2 = true_local_scores[1, 2],
    PC3 = true_local_scores[1, 3]
  )
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
pca_data <- cbind(pca_data_global |>
  dplyr::select(Design, Is_True_Map, PC1_global = PC1), 
  pca_data_local |> dplyr::select(PC1_local = PC1))

# Mahalanobis distance score (lower is better; true map should be near 0).
mahalanobis_safe <- function(x, center_point, tol = 1e-8) {
  x <- as.matrix(x)
  center_point <- as.numeric(center_point)
  cov_mat <- cov(x)

  # Robust inverse for singular/ill-conditioned covariance via eigenvalue flooring.
  eig <- eigen(cov_mat, symmetric = TRUE)
  eig_vals <- pmax(eig$values, tol)
  cov_inv <- eig$vectors %*% diag(1 / eig_vals) %*% t(eig$vectors)

  centered <- sweep(x, 2, center_point, FUN = "-")
  distances <- rowSums((centered %*% cov_inv) * centered)

  list(
    strategies = distances,
    center = center_point,
    cov_inv = cov_inv
  )
}

global_md <- mahalanobis_safe(metrics_global, center_point = rep(0, ncol(metrics_global)))
local_md <- mahalanobis_safe(metrics_local, center_point = rep(0, ncol(metrics_local)))

strategy_scores <- data.frame(
  Design = labels,
  score_global = sqrt(global_md$strategies),
  score_local = sqrt(local_md$strategies)
)

true_global_centered <- as.matrix(true_map_row_global) - matrix(global_md$center, nrow = 1)
true_local_centered <- as.matrix(true_map_row_local) - matrix(local_md$center, nrow = 1)
true_global_dist <- sqrt(rowSums((true_global_centered %*% global_md$cov_inv) * true_global_centered))
true_local_dist <- sqrt(rowSums((true_local_centered %*% local_md$cov_inv) * true_local_centered))

performance_scores <- rbind(
  strategy_scores,
  data.frame(
    Design = "True_Map",
    score_global = true_global_dist,
    score_local = true_local_dist,
    score_total = NA_real_
  )
)

# Normalize component scores so global and local/focal contribute equally.
normalize_component <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if ((rng[2] - rng[1]) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

performance_scores <- performance_scores |>
  mutate(
    score_global_norm = normalize_component(score_global),
    score_local_norm = normalize_component(score_local),
    score_total = sqrt(score_global_norm^2 + score_local_norm^2)
  ) |>
  arrange(score_total)

print(performance_scores)

performance_mean <- performance_scores %>%
  group_by(Design) %>%
  summarise(across(starts_with("score"), mean))

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

score_plot_data <- performance_scores |>
  mutate(Is_True_Map = ifelse(Design == "True_Map", "Yes", "No"))

p2 <- ggplot(score_plot_data, aes(x = score_global, y = score_local,
                                  color = Design, shape = Is_True_Map, size = Is_True_Map)) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = c("No" = 16, "Yes" = 8)) +
  scale_size_manual(values = c("No" = 3, "Yes" = 8)) +
  scale_color_manual(values = my_palette) +
  labs(x = "Global Mahalanobis distance", y = "Local/focal Mahalanobis distance") +
  theme_minimal() +
  guides(shape = "none", size = "none")

combined_plot <- p1 + p2 + 
  plot_annotation(
    title = "Multidimensional Performance of Sampling Designs",
    subtitle = "Left: PCA alignment between global and local/focal performance. Right: explicit distance-based error ranking."
  )
