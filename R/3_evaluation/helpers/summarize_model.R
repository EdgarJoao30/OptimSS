summarize_model <- function(fit, scenario) {
  # Extract fixed effects
  fixed_effects <- data.frame(
    Effect = rownames(fit$summary.fixed),
    Mean = fit$summary.fixed$mean,
    SD = fit$summary.fixed$sd,
    Type = "Fixed"
  )
  
  # Extract random effects
  random_effects <- data.frame(
    Effect = fit$summary.random$land_cover$ID,
    Mean = fit$summary.random$land_cover$mean,
    SD = fit$summary.random$land_cover$sd,
    Type = "Random"
  )
  
  # Extract hyperparameters
  hyperparameters <- data.frame(
    Effect = rownames(fit$summary.hyperpar),
    Mean = fit$summary.hyperpar$mean,
    SD = fit$summary.hyperpar$sd,
    Type = "Hyperparameter"
  )
  
  # timing 
  t <- sum(as.numeric(bru_timings(fit)$Elapsed))
  
  # Combine all effects into a single data frame
  combined_df <- rbind(
    cbind(fixed_effects, Scenario = scenario, time_s = t),
    cbind(random_effects, Scenario = scenario, time_s = t),
    cbind(hyperparameters, Scenario = scenario, time_s = t)
  )
  
  combined_df
}
