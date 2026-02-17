#################################################### Generate sampling scenarios ###############################################################
packs <- c("sf", "tidyverse")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
source('~/Documents/GitHub/OptimSS/R/2_sampling/helpers/sampling_scenarios.R')

roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')

sampling_designs <- generate_sampling_scenarios(roi,
                                   iterations = 100,
                                   sample_sizes = c(5, 10, 15, 25, 50),
                                   models = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'),
                                   crs = 32650)
samples <- sampling_designs$samples
grids <- sampling_designs$grids

st_write(samples, '~/Documents/GitHub/OptimSS/data/3_sampling/sampling_designs.geojson', append = FALSE)
# st_write(grids, paste0(wd, '/data/test/20250728_grids_sampling_scenarios_nobuffer.geojson'), append = F)
