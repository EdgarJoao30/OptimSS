library(INLA)
library(inlabru)
library(dplyr)
library(MASS)
### With INLA / inlabru

# Simulate 5000 points with land cover classes
set.seed(42)
n <- 5000
classes <- factor(sample(c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"), n, replace = TRUE))
levels(classes) <- c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built")

# Define contrasts with A_Primary as reference (baseline = 0)
true_beta <- c(
  "A_Primary" = 0,
  "B_Secondary" = 1.5,
  "C_Oil" = 1.0,
  "D_Plantation" = 0.5,
  "E_Built" = -0.5
)

# Simulate response
mu <- exp(true_beta[as.character(classes)])
y <- rpois(n, lambda = mu)

# Build data frame
dat <- data.frame(sim = y, class = classes)

# Fit Poisson model using inlabru
fit <- bru(sim ~ Intercept(1) + land_cover(class, model = "factor_contrast"),
           data = dat, family = "poisson")

# Predict expected value
pred <- predict(fit, dat, ~ exp(Intercept + land_cover))

# Compute pseudo-R²
poisson_loglik <- function(y, mu) {
  y * log(mu) - mu - lgamma(y + 1)
}

ll_model <- poisson_loglik(dat$sim, pred$mean)
ll_null <- poisson_loglik(dat$sim, rep(mean(dat$sim), n))

D_model <- -2 * sum(ll_model)
D_null <- -2 * sum(ll_null)

pseudo_r2 <- 1 - D_model / D_null
cat("Pseudo R²:", round(pseudo_r2, 3), "\n")



### With GLM 

set.seed(42)

# Simulate 5000 observations with land cover categories
n <- 5000
classes <- factor(sample(
  c("A_Primary", "B_Secondary", "C_Oil", "D_Plantation", "E_Built"),
  n,
  replace = TRUE
))

# Set true coefficients (on log scale), baseline = A_Primary (0)
# true_beta <- c(
#   "A_Primary" = 0,
#   "B_Secondary" = 1.5,
#   "C_Oil" = 1.0,
#   "D_Plantation" = 0.5,
#   "E_Built" = -0.5
# )

# true_beta <- c(
#   "A_Primary" = 0,
#   "B_Secondary" = 3,
#   "C_Oil" = 2,
#   "D_Plantation" = 1,
#   "E_Built" = -1
# )

true_beta <- c(
  "A_Primary" = 0,
  "B_Secondary" = 4.5,
  "C_Oil" = 3.0,
  "D_Plantation" = 1.5,
  "E_Built" = -1.5
)

# Linear predictor and simulate from Poisson
eta <- true_beta[as.character(classes)]
mu <- exp(eta)
# sim <- rpois(n, lambda = mu)
sim <- rnbinom(n = n, mu = mu, size = 100)

# Data frame
dat <- data.frame(sim = sim, class = classes)

# Fit GLM
# fit_glm <- glm(sim ~ class, data = dat, family = poisson(link = "log"))
fit_glm <- glm.nb(sim ~ class, data = dat)
summary(fit_glm)

# Predicted mean (mu) from the fitted model
mu_hat <- predict(fit_glm, type = "response")

# Null model (intercept only)
mu_null <- rep(mean(sim), n)

# Poisson log-likelihood function
poisson_loglik <- function(y, mu) {
  y * log(mu) - mu - lgamma(y + 1)
}

neg_bin_loglik <- function(y, mu, theta) {
  ll <- lgamma(y + theta) - lgamma(y + 1) - lgamma(theta) +
    theta * log(theta / (theta + mu)) +
    y * log(mu / (theta + mu))
  return(ll)
}

# Deviance calculation
# ll_model <- poisson_loglik(sim, mu_hat)
# ll_null <- poisson_loglik(sim, mu_null)
ll_model <- neg_bin_loglik(sim, mu_hat, fit_glm$theta)
ll_null <- neg_bin_loglik(sim, mu_null, fit_glm$theta)

D_model <- -2 * sum(ll_model)
D_null <- -2 * sum(ll_null)

pseudo_r2 <- 1 - D_model / D_null
cat("Pseudo R² (GLM):", round(pseudo_r2, 3), "\n")

plot(mu_hat, sim, pch = 20, col = scales::alpha("black", 0.3))
