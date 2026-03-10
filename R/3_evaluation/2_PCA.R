packs <- c("tidyverse")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

# load files
species <- commandArgs(trailingOnly = TRUE)[1]
species <- "aedes"
if (species == "anopheles") {
  range <- 3000
  rho <- 0.9
} else if (species == "aedes") {
  range <- 1000
  rho <- 0.7
} 

folder <- paste0("~/Documents/GitHub/OptimSS/data/4_results/results/", species, "/")
files <- list.files(folder, pattern = species, full.names = TRUE)
data <- lapply(files, read_csv) |> bind_rows()
true_data <- read.csv(paste0("~/Documents/GitHub/OptimSS/data/1_raw/", species, "_parameters.csv")) 

true_data <- true_data |>
  dplyr::select(ID, mean_abundance) |>
  spread(key = ID, value = mean_abundance)

colnames(true_data) <- c(
  "A_Primary_abundance",
  "B_Secondary_abundance",
  "C_Oil_abundance",
  "D_Plantation_abundance",
  "E_Built_abundance"
)

true_data <- true_data |>
  mutate(
    `size for the nbinomial observations (1/overdispersion)` = 10,
    `Range for field` = range,
    `Stdev for field` = 0.5,
    `GroupRho for field` = rho,
    global.global_abundance_rmse = 0,
    global.global_contig_rmse = 0,
    global.global_joinent_rmse = 0
  )


# remove rows with outliers
data <- data |> filter(is.na(focal.focal_abundance_rmse) | focal.focal_abundance_rmse < 15) |> 
  filter(is.na(`Range for field`) | `Range for field` < 5000)
labels <- data$scenario
sample_sizes <- data$sample_size
# Global metrics
metrics_global <- data |> dplyr::select(A_Primary_abundance, B_Secondary_abundance, C_Oil_abundance, D_Plantation_abundance, E_Built_abundance,
                                  `size for the nbinomial observations (1/overdispersion)`, `Range for field`, `Stdev for field`, `GroupRho for field`,
                                  global.global_abundance_rmse, global.global_contig_rmse, global.global_joinent_rmse)

metrics_global$A_Primary_abundance <- abs(metrics_global$A_Primary_abundance - true_data$A_Primary_abundance)
metrics_global$B_Secondary_abundance <- abs(metrics_global$B_Secondary_abundance - true_data$B_Secondary_abundance)
metrics_global$C_Oil_abundance <- abs(metrics_global$C_Oil_abundance - true_data$C_Oil_abundance)
metrics_global$D_Plantation_abundance <- abs(metrics_global$D_Plantation_abundance - true_data$D_Plantation_abundance)
metrics_global$E_Built_abundance <- abs(metrics_global$E_Built_abundance - true_data$E_Built_abundance)
metrics_global$`size for the nbinomial observations (1/overdispersion)` <- abs(metrics_global$`size for the nbinomial observations (1/overdispersion)` - true_data$`size for the nbinomial observations (1/overdispersion)`)
metrics_global$`Range for field` <- abs(metrics_global$`Range for field` - true_data$`Range for field`)
metrics_global$`Stdev for field` <- abs(metrics_global$`Stdev for field` - true_data$`Stdev for field`)
metrics_global$`GroupRho for field` <- abs(metrics_global$`GroupRho for field` - true_data$`GroupRho for field`)
# Local and focal metrics
metrics_local <- data |> dplyr::select(local.local_abundance_rmse, 
                                               focal.focal_abundance_rmse,
                                               focal.focal_contig_rmse,
                                               focal.focal_joinent_rmse, 
                                               focal.focal_classfreq_rmse)

# Missing estimates represent sampling limitations
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
  Sample_size = sample_sizes,
  PC1 = pca_res_global$x[, 1],
  PC2 = pca_res_global$x[, 2],
  PC3 = pca_res_global$x[, 3]
)
pca_data_global <- rbind(
  pca_data_global,
  data.frame(
    Design = "True_Map",
    Sample_size = 51,
    PC1 = true_global_scores[1, 1],
    PC2 = true_global_scores[1, 2],
    PC3 = true_global_scores[1, 3]
  )
)

pca_data_local <- data.frame(
  Design = labels,
  Sample_size = sample_sizes,
  PC1 = pca_res_local$x[, 1],
  PC2 = pca_res_local$x[, 2],
  PC3 = pca_res_local$x[, 3]
)
pca_data_local <- rbind(
  pca_data_local,
  data.frame(
    Design = "True_Map",
    Sample_size = 51,
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

pca_data <- cbind(pca_data_global |>
  dplyr::select(Design, Sample_size, Is_True_Map, PC1_global = PC1), 
  pca_data_local |> dplyr::select(PC1_local = PC1))



pca_data <- pca_data %>%
  mutate(space = case_when(
    Design == "True_Map" ~ "True_Map",
    Design %in% c("a", "d", "g") ~ "Lattice",
    Design %in% c("b", "e", "h") ~ "Stratified",
    Design %in% c("c", "f", "i") ~ "Random",
    TRUE ~ "Other"
  ),
  time = case_when(
    Design == "True_Map" ~ "True_Map",
    Design %in% c("a", "b", "c") ~ "Static",
    Design %in% c("d", "e", "f") ~ "Rotational",
    Design %in% c("g", "h", "i") ~ "Variable",
    TRUE ~ "Other"
  ))

write.csv(pca_data, paste0('~/Documents/GitHub/OptimSS/data/4_results/pca_data_', species, '.csv'), row.names = FALSE)
