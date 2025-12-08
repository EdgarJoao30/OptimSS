library(tidyverse)
library(sf)
library(INLA)
library(inlabru)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
mos <- read.csv(paste0(wd, '/data/mosquito_data_for_sim.csv')) %>% dplyr::select(-X.1)
mn <- seq(as.Date("2012/10/1"), by = "month", length.out = 123)
mn <- data.frame(mn, seq(1,123)); colnames(mn) <- c("date", "mon")
mos$date <- as.Date(paste0(mos$year,"/", mos$month, "/01"))
mos <- merge(mos, mn, by="date")
mos <- mos[c("date", "Class", "bite_rate", "X", "Y", "mon")]

mos <- st_as_sf(mos, coords = c('X', 'Y'), crs = 'EPSG: 4326') %>% st_transform(crs = 32650)
mos$Class <- plyr::revalue(mos$Class, c("Oil Palm"="Oil palm"))
mos$Class <- plyr::revalue(mos$Class, c("Oil palm"="C Oil palm", 
                                        "Primary" = "A Primary",
                                        "Secondary" = "B Secondary",
                                        "Plantation" = "D Plantation",
                                        "Kampung" = "E Built-up"))
#mos$Class <- relevel(factor(mos$Class), ref = "Primary")
# m1 <- MASS::glm.nb(formula = bite_rate ~ Class, data = mos)
# summary(m1)
mesh <- fm_mesh_2d(mos, max.edge = c(2500, 5000), cutoff = 1000)
# matern <- inla.spde2.matern(mesh, alpha = 2, constr = T)
matern <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.sigma = c(0.5, 0.01), prior.range = c(3000, 0.1))

# inlabru does not automatically remove the first level and absorb it into an intercept. 
# Instead, we can either use model="factor_full" without an intercept, or model="factor_contrast", 
# which does remove the first level.
# https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html

model <- bite_rate ~ Intercept(1) + 
  land_cover(Class, model = 'factor_contrast') + 
  field(geometry, model = matern, group = mon, control.group = list(model = 'ar1'))

fit <- bru(model, mos, family = "nbinomial",
           options = list(control.family = list(link = "log")
                          #control.compute = list(dic = TRUE, cpo = TRUE, config=T, dic = TRUE, waic = TRUE)
           ))

#save(fit, file = paste0(wd, "/data/models/estimate_parameters.RData"))

#load( paste0(wd, "/data/models/estimate_parameters.RData"))

sample <- inla.posterior.sample(10, result = fit)


flist <- vector("list", NROW(fit$summary.random$land_cover))
for (i in seq_along(flist)) flist[[i]] <- plot(fit, "land_cover", index = i)
multiplot(plotlist = flist, cols = 3)

spde.range <- spde.posterior(fit, "field", what = "range")
spde.logvar <- spde.posterior(fit, "", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot)
