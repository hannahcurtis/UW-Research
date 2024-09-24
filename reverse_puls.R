# Function to input a water level dataframe and get inflows and outflows

# make example reproducible
#set.seed(0)

# function to add a random value of the average noise to each water level measurement
add_noise <- function(elevation, ft_noise) {
  return(elevation + rnorm(length(elevation), 0, ft_noise))
}

# smoothing function (can change weights)
smooth_elevation <- function(elevation, weights=c(1/3,1/3,1/3)) {
  
  moving_avg_elevation <- rep(0,length(weights))
  for (i in 1:length(weights)) {
    j <- length(elevation) - (length(weights) - i)
    moving_avg_elevation <- moving_avg_elevation + weights[i]*elevation[i:j]
  }
  return(moving_avg_elevation)
}

# Transform the 1-minute data into 5-minute data
aggregate_data <- function(df, time_col, level_col, divisor) {
  # Calculate the number of complete groups of 3
  complete_rows <- (nrow(df) %/% divisor) * divisor
  
  # Filter out any incomplete group
  df <- df[1:complete_rows, ]
  
  df_xmin <- df %>%
    dplyr::mutate(group = rep(1:(n()/divisor), each = divisor)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      !!time_col := first(!!sym(time_col)),                # Use the first timestamp in each group
      !!level_col := mean(!!sym(level_col))               # Take the average of the three values
    ) %>%
    dplyr::select(-group)                                  # Drop the 'group' column
  
  return(df_xmin)
}

# Transform the 5-minute data into 15-minute data
aggregate_to_15min <- function(df, time_col, level_col) {
  # Calculate the number of complete groups of 3
  complete_rows <- (nrow(df) %/% 3) * 3
  
  # Filter out any incomplete group
  df <- df[1:complete_rows, ]
  
  df_15min <- df %>%
    dplyr::mutate(group = rep(1:(n()/3), each = 3)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      !!time_col := first(!!sym(time_col)),                # Use the first timestamp in each group
      !!level_col := mean(!!sym(level_col))               # Take the average of the three values
    ) %>%
    dplyr::select(-group)                                  # Drop the 'group' column
  
  return(df_15min)
}

# get storage difference by approximating pond as elliptical cylinder
get_storage_diff <- function(elevation, L_basin, W_basin) {
  # initialize storage_1 and storage_2
  storage_1 <- elevation
  storage_2 <- elevation
  # create a vector to store the differences in storage
  storage_diff <- numeric(length(elevation))
  for (i in 1:length(elevation)) {
    storage_1[i+1] <- storage_2[i] 
    storage_2[i+1] <- pi * elevation[i+1] * L_basin/2 * W_basin/2
    storage_diff[i] <- storage_2[i] - storage_1[i]
  }
  return(storage_diff)
}

# get storage difference by approximating pond as a truncated rectangular pyramid
get_storage_diff_pyr <- function(elevation, L_basin_long, W_basin_long, L_basin_short, W_basin_short) {
  # initialize storage_1 and storage_2
  storage_1 <- elevation
  storage_2 <- elevation
  # create a vector to store the differences in storage
  storage_diff <- numeric(length(elevation))
  for (i in 1:length(elevation)) {
    storage_1[i+1] <- storage_2[i] 
    storage_2[i+1] <- (((L_basin_long*W_basin_short) + (L_basin_short*W_basin_long) + (2*((L_basin_short*W_basin_short) + (L_basin_long*W_basin_long))))/6)*elevation[i+1]
    storage_diff[i] <- storage_2[i] - storage_1[i]
  }
  return(storage_diff)
}
  
# get average outflow for restrictive ponds
get_avg_outflow_res <- function(elevation, H_0, h_k, C_d, theta) {
  # initialize outflow_1 and outflow_2 with proper lengths
  outflow_1 <- numeric(length(elevation))
  outflow_2 <- numeric(length(elevation))
  # create a vector to store the average outflow
  avg_outflow <- numeric(length(elevation))
  # loop over the elevation vector
  for (i in 1:(length(elevation)-1)) {
    outflow_1[i+1] <- outflow_2[i]
    
    if (elevation[i] < H_0) {
      outflow_2[i+1] <- 0
    }
    else {
      outflow_2[i+1] <- (8/15)*C_d*sqrt(2*32.2)*tan(theta/2)*(((elevation[i]-H_0)+h_k)^(5/2))  
    }
    avg_outflow[i+1] <- (outflow_1[i+1] + outflow_2[i+1])/2
  }
  return(avg_outflow)
}

# get average outflow for non-restrictive ponds
get_avg_outflow_nonres <- function(elevation, H_0, L_w) {
  # initialize outflow_1 and outflow_2 with proper lengths
  outflow_1 <- numeric(length(elevation))
  outflow_2 <- numeric(length(elevation))
  # create a vector to store the average outflow
  avg_outflow <- numeric(length(elevation))
  for (i in 1:(length(elevation)-1)) {
    outflow_1[i+1] <- outflow_2[i]
    
    if (elevation[i] < H_0) {
      outflow_2[i+1] <- 0
    }
    else {
      outflow_2[i+1] <- 3.0888*L_w*((elevation[i]-H_0)^(3/2))
    }
  avg_outflow[i+1] <- (outflow_1[i+1] + outflow_2[i+1])/2  
  }
  return(avg_outflow)
}

# get average outflow for EPA-SWMM runs
get_avg_outflow_swmm <- function(elevation, H_0, C_w, SS) {
  # initialize outflow_1 and outflow_2 with proper lengths
  outflow_1 <- numeric(length(elevation))
  outflow_2 <- numeric(length(elevation))
  # create a vector to store the average outflow
  avg_outflow <- numeric(length(elevation))
  # loop over the elevation vector
  for (i in 1:(length(elevation)-1)) {
    outflow_1[i+1] <- outflow_2[i]
    
    if (elevation[i] < H_0) {
      outflow_2[i+1] <- 0
    }
    else {
      outflow_2[i+1] <- (C_w*SS)*((elevation[i])^(4.1/2))  
    }
    avg_outflow[i+1] <- (outflow_1[i+1] + outflow_2[i+1])/2
  }
  return(avg_outflow)
}

# get average inflow 
get_avg_inflow <- function(storage_diff, avg_outflow, time_delta) {
  return((storage_diff/time_delta)+avg_outflow)
}

# calculates the inflow and outflow for restrictive ponds
get_inflow_outflow_res <- function(elevation, ft_bias, ft_noise, L_basin, W_basin, H_0, h_k, C_d, theta, time_delta) {
  # add noise to each elevation measurement
  noise_elevation <- add_noise(elevation, ft_noise)
  # add bias to each elevation measurement
  bias_elevation <- noise_elevation + ft_bias
  # calculate storage diff
  storage_diff <- get_storage_diff(bias_elevation, L_basin, W_basin)
  # calculate average outflow
  avg_outflow <- get_avg_outflow_res(bias_elevation, H_0, h_k, C_d, theta)
  # calculate average inflow
  avg_inflow <- get_avg_inflow(storage_diff, avg_outflow, time_delta)
  # create dataframe with inflow and outflow
  result_df <- data.frame(Inflow.cfs = avg_inflow, Outflow.cfs = avg_outflow)
  # return the dataframe
  return(result_df)
}

# calculates the inflow and outflow for non-restrictive ponds
get_inflow_outflow_nonres <- function(elevation, ft_bias, ft_noise, L_basin, W_basin, H_0, L_w, time_delta) {
  # add noise to each elevation measurement
  noise_elevation <- add_noise(elevation, ft_noise)
  # add bias to each elevation measurement
  bias_elevation <- noise_elevation + ft_bias
  # calculate storage diff
  storage_diff <- get_storage_diff(bias_elevation, L_basin, W_basin)
  # calculate average outflow
  avg_outflow <- get_avg_outflow_nonres(bias_elevation, H_0, L_w)
  # calculate average inflow
  avg_inflow <- get_avg_inflow(storage_diff, avg_outflow, time_delta)
  # create dataframe with inflow and outflow
  result_df <- data.frame(Inflow.cfs = avg_inflow, Outflow.cfs = avg_outflow)
  # return the dataframe
  return(result_df)
}

get_inflow_outflow_swmm <- function(elevation, ft_bias, ft_noise, L_basin, W_basin, H_0, C_w, SS, time_delta) {
  # add noise to each elevation measurement
  noise_elevation <- add_noise(elevation, ft_noise)
  # add bias to each elevation measurement
  bias_elevation <- noise_elevation + ft_bias
  # calculate storage diff
  storage_diff <- get_storage_diff(bias_elevation, L_basin, W_basin)
  # calculate average outflow
  avg_outflow <- get_avg_outflow_swmm(bias_elevation, H_0, C_w, SS)
  # calculate average inflow
  avg_inflow <- get_avg_inflow(storage_diff, avg_outflow, time_delta)
  # create dataframe with inflow and outflow
  result_df <- data.frame(Inflow.cfs = avg_inflow, Outflow.cfs = avg_outflow)
  # return the dataframe
  return(result_df) 
}

calc_magnitude_difference <- function(original, new) {
  return(original - new)
}

calc_percent_difference <- function(original, new) {
  return(((original - new)/original)*100)
}


