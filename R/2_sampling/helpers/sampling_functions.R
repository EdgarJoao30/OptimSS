#################################################### Functions for sampling every month ###############################################################
library(sf)
library(sp)
library(raster)
library(tidyverse)

# wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
# roi <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |>
#   st_transform(crs = 32650)


############## 
############## Function to sample randomly 
############## 

random_sample <- function(roi, sample_size = 15, points_per_sample = 3, fixed = TRUE) {
  # If fixed, sample grids once for the entire year
  if (fixed) {
    sampled_indices <- sample(1:nrow(roi), sample_size)
    sampled_grids <- roi[sampled_indices, ]
  }
  
  # Initialize lists to store samples and grids for all months
  all_samples <- vector("list", 12)
  all_grids <- vector("list", 12)
  
  # Loop through each month
  for (month in 1:12) {
    # If not fixed, sample grids for each month
    if (!fixed) {
      sampled_indices <- sample(1:nrow(roi), sample_size)
      sampled_grids <- roi[sampled_indices, ]
    }
    
    # Generate random points within each sampled grid
    month_samples <- lapply(1:nrow(sampled_grids), function(i) {
      pts <- spsample(as(sampled_grids[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- i
      return(pts_df)
    })
    
    # Combine all samples for the month
    all_samples[[month]] <- do.call(rbind, month_samples)
    
    # Add the month column to the grid data
    sampled_grids$month <- month
    all_grids[[month]] <- sampled_grids
  }
  
  # Combine all months into single data frames
  combined_samples <- do.call(rbind, all_samples)
  combined_grids <- do.call(rbind, all_grids)
  
  return(list(combined_samples, combined_grids))
}

# test <- uniform_sample(roi, fixed = FALSE)
# 
# test_samples <- st_as_sf(test[[1]], coords = c('x', 'y'), crs = 32650)
# test_grids <- test[[2]]
# 
# ggplot() +
#   #geom_sf(data = roi) +
#   geom_sf(data = test_grids %>% filter(month == 12)) +
#   geom_sf(data = test_samples%>% filter(month == 12), aes(color = month))
# 


############## 
############## Function to make a stratified sample using the categories 
############## 

stratified_sample <- function(roi, sample_size = 15, points_per_sample = 3, fixed = TRUE) {
  # Initialize lists to store samples and grids for all months
  all_grids <- vector("list", 12)
  all_samples <- vector("list", 12)
  
  # If fixed, sample grids once for the entire year
  if (fixed) {
    categories <- unique(roi$cat)
    samples_per_category <- ceiling(sample_size / length(categories))
    
    fixed_grids <- lapply(categories, function(cat) {
      cat_shp <- roi[roi$cat == cat, ]
      sampled_indices <- sample(1:nrow(cat_shp), min(samples_per_category, nrow(cat_shp)))
      cat_shp[sampled_indices, ]
    })
    fixed_grids <- do.call(rbind, fixed_grids)
  }
  
  # Loop through each month
  for (month in 1:12) {
    if (!fixed) {
      # Sample grids for the current month
      categories <- unique(roi$cat)
      samples_per_category <- ceiling(sample_size / length(categories))
      
      month_grids <- lapply(categories, function(cat) {
        cat_shp <- roi[roi$cat == cat, ]
        sampled_indices <- sample(1:nrow(cat_shp), min(samples_per_category, nrow(cat_shp)))
        cat_shp[sampled_indices, ]
      })
      month_grids <- do.call(rbind, month_grids)
    } else {
      # Use the fixed grids for all months
      month_grids <- fixed_grids
    }
    
    # Generate random points within each sampled grid
    month_samples <- lapply(1:nrow(month_grids), function(i) {
      pts <- spsample(as(month_grids[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- paste(month_grids$cat[i], i, sep = "_")
      return(pts_df)
    })
    
    # Combine all samples for the month
    all_samples[[month]] <- do.call(rbind, month_samples)
    
    # Add the month column to the grid data
    month_grids$month <- month
    all_grids[[month]] <- month_grids
  }
  
  # Combine all months into single data frames
  combined_samples <- do.call(rbind, all_samples)
  combined_grids <- do.call(rbind, all_grids)
  
  return(list(combined_samples, combined_grids))
}

############## 
############## Function to get samples equally spaced in the overall geometry of all the polygons together for each month
############## 

uniform_sample <- function(roi, sample_size = 15, points_per_sample = 3, fixed = TRUE) {
  # Initialize lists to store samples and grids for all months
  all_samples <- vector("list", 12)
  all_grids <- vector("list", 12)
  
  # If fixed, generate equally spaced sample locations once for the entire year
  if (fixed) {
    dissolved <- roi %>% st_union()
    
    repeat {
      samples_loc <- st_sample(dissolved, size = sample_size, type = 'regular')
      if (length(samples_loc) == sample_size) {
        break
      }
    }
    
    fixed_joined <- st_join(roi, st_as_sf(samples_loc), join = st_intersects, left = FALSE)
  }
  
  # Loop through each month
  for (month in 1:12) {
    if (!fixed) {
      # Generate equally spaced sample locations for the current month
      dissolved <- roi %>% st_union()
      
      repeat {
        samples_loc <- st_sample(dissolved, size = sample_size, type = 'regular')
        if (length(samples_loc) == sample_size) {
          break
        }
      }
      
      month_joined <- st_join(roi, st_as_sf(samples_loc), join = st_intersects, left = FALSE)
    } else {
      # Use the fixed sample locations for all months
      month_joined <- fixed_joined
    }
    
    # Generate random points within each sampled grid
    month_samples <- lapply(1:nrow(month_joined), function(i) {
      pts <- spsample(as(month_joined[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- i
      return(pts_df)
    })
    
    # Combine all samples for the month
    all_samples[[month]] <- do.call(rbind, month_samples)
    
    # Add the month column to the grid data
    month_joined$month <- month
    all_grids[[month]] <- month_joined
  }
  
  # Combine all months into single data frames
  combined_samples <- do.call(rbind, all_samples)
  combined_grids <- do.call(rbind, all_grids)
  
  return(list(combined_samples, combined_grids))
}
