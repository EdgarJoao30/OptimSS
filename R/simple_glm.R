pacman::p_load(INLA, inlabru, sf, fmesher, terra, tidyverse)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif'))

points <- as.data.frame(landcover, xy=T)
points$class <- factor(points$Landcover_AllClass)
points$class <- plyr::revalue(points$class, c("0"="C_Oil", "1"="B_Secondary", "2"="A_Primary", "3"="D_Plantation", "4"="E_Built"))
points$class <- factor(points$class, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
points <- st_as_sf(points, coords = c('x', 'y'), crs = 32650)

### SIMULATION
generate_nbinomial <- function(x) {
  rnbinom(mu = x, n = 1, size = 10)
}
beta <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367) # Primary, Secondary, Oil Palm, Plantation, Built-up
beta0 <- beta[1]
beta[1] <- 0
months <- 1
ccov <- factor(replicate(months, points$class))
mu <- beta0 + beta[unclass(ccov)]
points$mu <- exp(mu)
seed <- 1234
set.seed(seed)
nbinomial_sample <- sapply(points$mu, generate_nbinomial)
points$sim <- nbinomial_sample
names(points)[1] <- 'landcover'

### MODEL FIT
sample_points <- points %>% slice_sample(n = 1000)
model <- sim ~ Intercept(1)  + 
  land_cover(class, model = 'factor_contrast') 

fit <- bru(
  model,
  sample_points,
  family = "nbinomial",
  options = list(
    control.family = list(link = "log"),
    control.inla = list(int.strategy = "eb")
    # control.fixed = list(
    #   mean = 0,
    #   prec = 1e-6  # Very large variance = flat prior
    # )
  )
)

summary(fit)
fit$summary.fixed
fit$summary.random$land_cover

# PREDICT
pred <- 
  predict(fit, points,
          ~ {
            mu <- exp(Intercept + land_cover)
            pred <- rnbinom(n = nrow(points), size = fit$summary.hyperpar$mean[1], mu = mu)
            
            list(
              mu = mu,
              pred = pred
            )
          },
          n.samples = 1000
  )

### EVALUATE
neg_bin_loglik <- function(y, mu, theta) {
  ll <- lgamma(y + theta) - lgamma(y + 1) - lgamma(theta) +
    theta * log(theta / (theta + mu)) +
    y * log(mu / (theta + mu))
  return(ll)
}
y <- points$sim
mu_model <- pred$pred$mean
mu_null <- rep(mean(y), length(y))
k <- fit$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
ll_model <- neg_bin_loglik(y, mu_model, k)
ll_null <- neg_bin_loglik(y, mu_null, k)

D_model <- -2 * sum(ll_model)
D_null <- -2 * sum(ll_null)

pseudo_r2 <- 1 - D_model / D_null

pred$pred %>% 
ggplot() +
  geom_point(
    aes(
      x = sim,
      y = mean,
      fill = as.factor(landcover)
    ), pch = 21, alpha = .5
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("#f6e8c3", "#80cdc1", "#018571", "#d8b365", "gray"),
    labels = c(
      "0" = "Oil palm.",
      "1" = "Secondary forest.",
      "2" = "Primary forest.",
      "3" = "Other plantations.",
      "4" = "Built-up."
    ),
    guide = guide_legend(position = "bottom", ncol = 2)
  ) +
  theme_minimal(base_size = 12, base_family = 'Times New Roman') +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(size = rel(1.5), face = "italic")
  ) +
  labs(x = 'Simulation', y = 'Prediction', fill = 'Land cover')
