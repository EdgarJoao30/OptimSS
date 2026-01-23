### Load and install required packages
packs <- c("tidyverse", "INLA", "inlabru", "sf")
success <- suppressWarnings(sapply(packs, require, character.only = TRUE))
install.packages(names(success)[!success])
sapply(names(success)[!success], require, character.only = TRUE)
### Load data
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'

an_raw <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/anopheles_raw.geojson')
ae_raw <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/aedes_raw.geojson')
roi <- st_read('~/Documents/GitHub/OptimSS/data/1_raw/ROI_32650.geojson')
mos <- ae_raw
### Hyperparameters
range <- 1000
rho <- 0.9
### model fit
mesh <- fm_mesh_2d(roi, max.edge = c(2500, 5000), cutoff = 1000)
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, 0.01), prior.range = c(range, 0.1))
model <- bite_rate ~ Intercept(1) + 
  land_cover(Class, model = 'factor_contrast') + 
  field(geometry, model = matern, group = mon, 
        control.group = list(model = 'ar1', 
                             hyper = list(theta=list(prior='normal', param=c(rho, 0.001))) ))

fit <- bru(model, mos, family = "nbinomial",
           options = list(control.family = list(link = "log"),
                          control.inla = list(int.strategy = "eb")
           ))
summary(fit)
fit$summary.random$land_cover

# save params in a csv ready to be read by next script