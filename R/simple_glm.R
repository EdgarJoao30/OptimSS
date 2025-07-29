pacman::p_load(INLA, inlabru, sf, fmesher, terra, tidyverse)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif'))

points <- as.data.frame(landcover, xy=T)
points$class <- factor(points$Landcover_AllClass)
points$class <- plyr::revalue(points$class, c("0"="C_Oil", "1"="B_Secondary", "2"="A_Primary", "3"="D_Plantation", "4"="E_Built"))
points$class <- factor(points$class, levels = c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"))
points <- st_as_sf(points, coords = c('x', 'y'), crs = 32650)
names(points)[1] <- 'landcover'

# mos <- read.csv(paste0(wd, '/data/mosquito_data_for_sim.csv')) %>% dplyr::select(-X.1)
# mn <- seq(as.Date("2012/10/1"), by = "month", length.out = 123)
# mn <- data.frame(mn, seq(1,123)); colnames(mn) <- c("date", "mon")
# mos$date <- as.Date(paste0(mos$year,"/", mos$month, "/01"))
# mos <- merge(mos, mn, by="date")
# mos <- mos[c("date", "Class", "bite_rate", "X", "Y", "mon")]
# 
# mos <- st_as_sf(mos, coords = c('X', 'Y'), crs = 'EPSG: 4326') %>% st_transform(crs = 32650)
# mos$class <- plyr::revalue(mos$Class, c("Oil Palm"="Oil palm"))
# mos$class <- plyr::revalue(mos$Class, c("Oil palm"="C_Oil", 
#                                         "Primary" = "A_Primary",
#                                         "Secondary" = "B_Secondary",
#                                         "Plantation" = "D_Plantation",
#                                         "Kampung" = "E_Built"))


### SIMULATION
generate_nbinomial <- function(x) {
  rnbinom(mu = x, n = 1, size = 10)
}

generate_poisson <- function(x) {
  rpois(n = 1, lambda = x)
}
beta <- c(0.7351552, 2.588883, 1.775618, 1.612111, 1.156367) # Primary, Secondary, Oil Palm, Plantation, Built-up
beta0 <- beta[1]
beta[1] <- 0

# Set the intercept to 0
beta0 <- 0
# Set A_Primary as reference (0), and use contrasts relative to it
# So beta = 0 for A_Primary, and other values are relative differences
beta_contrasts <- c(
  "A_Primary" = 0,
  "B_Secondary" = 2.588883 - 0.7351552,
  "C_Oil" = 1.775618 - 0.7351552,
  "D_Plantation" = 1.612111 - 0.7351552,
  "E_Built" = 1.156367 - 0.7351552
)

points$mu <- exp(beta0 + beta_contrasts[as.character(points$class)])
points$sim <- sapply(points$mu, generate_poisson)

months <- 1
ccov <- factor(replicate(months, points$class))
mu <- beta0 + beta[unclass(ccov)]
points$mu <- exp(mu)
seed <- 1234
set.seed(seed)
dist_sample <- sapply(points$mu, generate_poisson)
points$sim <- dist_sample



# model <- bite_rate ~ Intercept(1)  + 
#   land_cover(class, model = 'factor_contrast') 
# 
# fit <- bru(
#   model,
#   mos,
#   family = "nbinomial",
#   options = list(
#     control.family = list(link = "log"),
#     control.inla = list(int.strategy = "eb")
#     # control.fixed = list(
#     #   mean = 0,
#     #   prec = 1e-6  # Very large variance = flat prior
#     # )
#   )
# )
# 
# summary(fit)
# fit$summary.fixed
# fit$summary.random$land_cover
# 
# sim <- generate(fit, points,
#          ~ {
#            mu <- exp(Intercept + land_cover)
#            sim <- rnbinom(n = nrow(points), size = fit$summary.hyperpar$mean[1], mu = mu)
#            },
#          n.samples = 1
# )
# 
# points$sim <- sim[, 1]


### MODEL FIT
sample_points <- points %>% slice_sample(n = 1000)
model <- sim ~ Intercept(1)  + 
  land_cover(class, model = 'factor_contrast') 

fit <- bru(
  model,
  sample_points,
  family = "poisson",
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
            # pred <- rnbinom(n = nrow(points), size = fit$summary.hyperpar$mean[1], mu = mu)
            pred <- sapply(mu, generate_poisson)
            
            list(
              mu = mu,
              pred = pred
            )
          },
          n.samples = 1000
  )

# Use only sample_points for evaluation
pred_sample <- predict(fit, sample_points,
                       ~ {
                         mu <- exp(Intercept + land_cover)
                         # pred <- rnbinom(n = nrow(points), size = fit$summary.hyperpar$mean[1], mu = mu)
                         pred <- sapply(mu, generate_poisson)
                         
                         list(
                           mu = mu,
                           pred = pred
                         )
                       },
                       n.samples = 1000
)

### EVALUATE

poisson_loglik <- function(y, mu) {
  ll <- y * log(mu) - mu - lgamma(y + 1)
  return(ll)
}

neg_bin_loglik <- function(y, mu, theta) {
  ll <- lgamma(y + theta) - lgamma(y + 1) - lgamma(theta) +
    theta * log(theta / (theta + mu)) +
    y * log(mu / (theta + mu))
  return(ll)
}
y <- points$sim
mu_model <- pred$pred$mean
mu_null <- rep(mean(y), length(y))

y_sample <- sample_points$sim
mu_model_sample <- pred_sample$pred$mean
mu_null_sample <- rep(mean(y_sample), length(y_sample))
# k <- fit$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
ll_model <- poisson_loglik(y, mu_model)
ll_null <- poisson_loglik(y, mu_null)
ll_model_sample <- poisson_loglik(y_sample, mu_model_sample)
ll_null_sample <- poisson_loglik(y_sample, mu_null_sample)
D_model <- -2 * sum(ll_model)
D_null <- -2 * sum(ll_null)
D_model_sample <- -2 * sum(ll_model_sample)
D_null_sample <- -2 * sum(ll_null_sample)
pseudo_r2 <- 1 - D_model / D_null
pseudo_r2_sample <- 1 - D_model_sample / D_null_sample


pred_sample$pred %>% 
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
