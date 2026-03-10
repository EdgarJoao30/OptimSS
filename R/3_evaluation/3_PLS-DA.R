# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("mixOmics")

packs <- c("tidyverse", "mixOmics")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)

# load files
species <- commandArgs(trailingOnly = TRUE)[1]
species <- "anopheles"
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
labels <- as.factor(data$scenario)
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


plsda_global <- plsda(metrics_global, labels, ncomp = 2)