library(tidyverse)
library(sf)
library(INLA)
library(inlabru)
library(lubridate)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
mos <- read.csv(paste0(wd, '/data/aedes_xy.csv')) %>% 
  mutate(date = dmy(date)) %>% 
  st_as_sf(coords = c('X', 'Y'), crs = 'EPSG: 4326') %>% 
  st_transform(crs = 32650)
mos$class <- plyr::revalue(mos$class, c("oil"="C Oil palm", 
                                        "primary" = "A Primary",
                                        "secondary" = "B Secondary",
                                        "plantation" = "D Plantation",
                                        "built" = "E Built-up"))
roi <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_transform(crs = 32650) #%>% st_buffer(-10)

ggplot()+
  geom_sf(data = roi) +
  geom_sf(data = mos)

#mos$Class <- relevel(factor(mos$Class), ref = "Primary")
m1 <- MASS::glm.nb(formula = bite_rate ~ class, data = mos)
summary(m1)
mesh <- fm_mesh_2d(roi, max.edge = c(2500, 5000), cutoff = 1000)
# matern <- inla.spde2.matern(mesh, alpha = 2, constr = T)
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, 0.01), prior.range = c(1000, 0.1))

# inlabru does not automatically remove the first level and absorb it into an intercept. 
# Instead, we can either use model="factor_full" without an intercept, or model="factor_contrast", 
# which does remove the first level.
# https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html

model <- bite_rate ~ Intercept(1) + 
  land_cover(class, model = 'factor_contrast') + 
  field(geometry, model = matern, group = month, control.group = list(model = 'ar1'))

fit <- bru(model, mos, family = "nbinomial",
           options = list(control.family = list(link = "log"),
                          control.inla = list(int.strategy = "eb")
           ))

#save(fit, file = paste0(wd, "/data/models/ae_estimate_parameters.RData"))

#load( paste0(wd, "/data/models/ae_estimate_parameters.RData"))

flist <- vector("list", NROW(fit$summary.random$land_cover))
for (i in seq_along(flist)) flist[[i]] <- plot(fit, "land_cover", index = i)
multiplot(plotlist = flist, cols = 3)

spde.range <- spde.posterior(fit, "field", what = "range")
spde.logvar <- spde.posterior(fit, "field", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot)
