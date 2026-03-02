summarize_model <- function(fit, scenario, n_levels = NULL) {
  # Extract fixed effects
  fixed_effects <- data.frame(
    effect = rownames(fit$summary.fixed),
    mean = fit$summary.fixed$mean,
    sd = fit$summary.fixed$sd,
    type = "Fixed"
  )
  
  # Extract random effects
  if (!is.null(n_levels) && n_levels > 1) {
    random_effects <- data.frame(
      effect = fit$summary.random$land_cover$ID,
      mean = fit$summary.random$land_cover$mean,
      sd = fit$summary.random$land_cover$sd,
      type = "Random"
    )
  }
  
  # Extract hyperparameters
  hyperparameters <- data.frame(
    effect = rownames(fit$summary.hyperpar),
    mean = fit$summary.hyperpar$mean,
    sd = fit$summary.hyperpar$sd,
    type = "Hyperparameter"
  )
  
  # timing 
  t <- sum(as.numeric(bru_timings(fit)$Elapsed))
  
  # Combine all effects into a single data frame
  if (!is.null(n_levels) && n_levels > 1) {
    combined_df <- rbind(
      cbind(fixed_effects, scenario = scenario, time_s = t),
      cbind(random_effects, scenario = scenario, time_s = t),
      cbind(hyperparameters, scenario = scenario, time_s = t)
    )
  } else {
    combined_df <- rbind(
      cbind(fixed_effects, scenario = scenario, time_s = t),
      cbind(hyperparameters, scenario = scenario, time_s = t)
    )
  }
  combined_df
}
