suppressMessages(library(tidyverse))
suppressMessages(library(terra))
suppressMessages(library(sf))

ssdwd <- '/Volumes/Anopheles/PhD'
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
files_dir <- paste0(ssdwd, '/post_samples')


ss <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
element_i <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
#ss <- 50

element <- letters[element_i]

files <- list.files(files_dir, full.names = TRUE, pattern = paste0("_s", ss, "_"))

sampled_data_list <- list()



for (file in files) {
    print(file)
    load(file)
    print(element)
    if (!is.null(sampled_data_list[[element]])) {
        sampled_data <- matrix(nrow = 72336, ncol = 100)
        for (i in 1:72336) {
            sampled_data[i, ] <- sample(samp_list[[element]][i, ], size = 100, replace = FALSE)
        }
        sampled_data_list[[element]] <- cbind(sampled_data_list[[element]], sampled_data)
    } else {
        sampled_data <- matrix(nrow = 72336, ncol = 100)
        for (i in 1:72336) {
            sampled_data[i, ] <- sample(samp_list[[element]][i, ], size = 100, replace = FALSE)
        }
        sampled_data_list[[element]] <- sampled_data
    }
}

save(sampled_data_list, file = paste0(files_dir, "/lite/", "post_", element, '_', ss, ".RData"))
