source("tools/load_libraries.R")

# Global variables
list_mds <<- list()
n_recursive_calls <<- 1
list_sample_points <<- list()
M_align <<- data.frame()
n_obs_main_matrix <<- NA
n_sampling_points <<- NA


# This function decides whether it is possible to compute distance matrix
is_possible_to_calculate_distance_matrix <- function(
  x,
  metric,
  timeout
  
){
  
  plan(multisession, gc = TRUE)
  
  
  # Calculate distance
  res_distance = NULL
  res_distance %<-% {
    daisy(
      x = x,
      metric = metric
    )
  }
  
  is_able_to_resolve = NULL
  is_able_to_resolve <- withTimeout(
    {
      resolved(res_distance)
    },
    timeout = timeout,
    onTimeout = "silent"
  )
  
  is_distance_computed = is.null(is_able_to_resolve) == FALSE && 
    is.na(is_able_to_resolve) == FALSE &&
    is_able_to_resolve == TRUE
  
  list_to_return = list(
    is_distance_computed = is_distance_computed,
    distance_matrix = res_distance
  )
  
  gc()
  
  return(list_to_return)
}

# This function decides whether it is possible to compute MDS
is_possible_to_run_mds <- function(
  distance_matrix,
  number_coordinates,
  timeout
){
  plan(multisession, gc = TRUE)
  
  # Calculating MDS
  mds_classical %<-% {
    stats::cmdscale(
      d = distance_matrix, 
      k = number_coordinates
    )
  }
    
  is_able_to_resolve = NULL
  is_able_to_resolve <- withTimeout(
    {
      resolved(mds_classical)
    },
    timeout = timeout,
    onTimeout = "silent"
  )
    
  is_mds_computed = is.null(is_able_to_resolve) == FALSE && 
    is.na(is_able_to_resolve) == FALSE &&
    is_able_to_resolve == TRUE
  
  
  gc()
  
  return(is_mds_computed)
    
}
  
  

# This function decides whether it is possible to compute distance and MDS
is_possible_to_calculate_components_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout
){
  
  # Check if we can compute distance matrix
  result_distance_computation = is_possible_to_calculate_distance_matrix(
    x = x,
    metric = metric,
    timeout = timeout
  )
  
  
  able_to_compute_distance = result_distance_computation$is_distance_computed
  
  
  # Check if we can compute MDS
  if( able_to_compute_distance == TRUE ){
    able_to_compute_components_mds = is_possible_to_run_mds(
      distance_matrix = result_distance_computation$distance_matrix,
      number_coordinates = number_coordinates,
      timeout = timeout
    )
    
  }
  
  gc()
  return(able_to_compute_components_mds)
}


recursive_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout
){
  
  # Store the rowns of the main data frame
  if( is.na(n_obs_main_matrix) == TRUE ){
    n_obs_main_matrix <<- nrow(x) 
    message(paste0("Rows of the main matrix: ", n_obs_main_matrix))
  } 
  
  
  # Amplify the number of coordinates by a factor
  amplification_factor = 3
  number_coordinates_amplified = amplification_factor * number_coordinates
  
  # Check if it is possible to run MDS
  is_possible_mds = is_possible_to_calculate_components_mds(
    x = x,
    number_coordinates = number_coordinates_amplified,
    metric = metric,
    timeout = timeout
  )
  
  # It is possible to run the MDS, the do it
  if(
    is.na(is_possible_mds) == FALSE && 
    is.null(is_possible_mds) == FALSE && 
    is_possible_mds == TRUE
  ){
    message(paste0("Possible to run MDS with ", nrow(x), " rows"))
    message(paste0("Value of n_recursive_calls: ", n_recursive_calls))
    
    
    # calculate the sampling points to compute M_align
    if( is.na(n_sampling_points) == TRUE){
      p = ceiling(n_obs_main_matrix/nrow(x))
      n_sampling_points <<- floor(nrow(x)/p)
      message(paste0("Number of sampling points: ", n_sampling_points))
    }
    
    # Calculate distance
    distance_matrix = daisy(
      x = x,
      metric = metric
    )
    
    # Applying MDS to the distance matrix
    classical_mds = stats::cmdscale(
      d = distance_matrix, 
      k = number_coordinates
    ) 
    
    row.names(classical_mds) = row.names(x)
    
    output_mds = classical_mds
    
    # Save the output
    list_mds[[n_recursive_calls]] <<- output_mds
    
    # Take s points from the matrix randomly
    samled_points = sample(
      x = row.names(x), 
      size = n_sampling_points, 
      replace = FALSE
    )
    
    list_sample_points[[n_recursive_calls]] <<- samled_points
    
    # Append the points no M_align
    M_align <<- rbind(M_align, x[samled_points, ])

    
    # Increment for the next MDS
    n_recursive_calls <<- n_recursive_calls + 1
  
  }else{
    message(paste0("Unable to run MDS for ", nrow(x), " observations"))
    # Divide the data into two parts and run call recursively
    n_division = 2
    group_distribution = sample(x = n_division, size = nrow(x), replace = TRUE) 
    for(i_group in 1:n_division){
      # Splitting x into small matrices to call MDS
      ind = which(group_distribution == i_group)
      x_divide = x[ind, ]
      
      # Call to the fast MDS
      recursive_mds(
        x = x_divide,
        number_coordinates = number_coordinates,
        metric = metric,
        timeout =  timeout
      )
    }
  }
}

fast_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout = 60

){
  #Calling the recursive process
  recursive_mds(
    x = x,
    number_coordinates = number_coordinates,
    metric = metric,
    timeout = timeout
  )
  message("Finished the recursive MDS")
  
  # Once it is done, we stitch the solutions so that all the points have the
  # same system of coordinates. 
  
  # Perform MDS on M_align
  distance_M_align = daisy(
    x = M_align,
    metric = metric
  )
  
  mds_M_align = stats::cmdscale(
    d = distance_M_align, 
    k = number_coordinates
  )
  
  # We stitch the solutions
  total_groups = length(list_mds)
  
  for(i_group in 1:total_groups){
    # Filter the s points from D_i and M_align
    rows_to_select = list_sample_points[[i_group]]
    d_mds_i = list_mds[[i_group]]
    d_mds_i_filter = d_mds_i[rows_to_select,] 
    M_align_filter = mds_M_align[rows_to_select, ]
    
    # Solving procruster problem
    procrustes_result = smacof::Procrustes(
      X = M_align_filter, 
      Y = d_mds_i_filter
    )
    
    rotation_matrix = procrustes_result$rotation
    dilation = procrustes_result$dilation
    translation = procrustes_result$translation
    
    
    # Transforming the data for the k-th group  
    tranformation_di = dilation * d_mds_i %*% rotation_matrix + translation
    
    if(i_group == 1){
      mds_aligned = tranformation_di

    }else{
      mds_aligned = rbind(
        mds_aligned,
        tranformation_di
      )
    }
    
  }
  
  
  # Reordering
  reording_permutation = match(row.names(x), row.names(mds_aligned))
  mds_aligned = mds_aligned[reording_permutation, ]
  
  
  return(mds_aligned)
  
}


###########################################################################
# Checks for fast_mds function
n_obs = 10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)

results_fast_mds = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = "euclidean",
  timeout = 1
)



head(results_fast_mds, 10)
dim(M_align)
length(list_mds)

distance_x = daisy(
  x = x,
  metric = "euclidean"
)


mds_x = stats::cmdscale(
  d = distance_x, 
  k = 2
)



rotation_matrix = smacof::Procrustes(
  X = mds_x,
  Y = results_fast_mds
)


df_classical = as.data.frame(mds_x)
df_classical$type = "classical"
df_classical$label = row.names(x)


df_fast = as.data.frame(rotation_matrix$Yhat)
df_fast$type = 'divide_conquer'
df_fast$label = row.names(x)

df_all = rbind(
  df_classical,
  df_fast
) %>% 
  arrange(
    as.numeric(label)
  )

# Plot
df_all[1:20,] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Euclidean distance
distance_classical_divide = diag(
  rdist(
    mds_classical,
    procrustes_result$Yhat
  )
)

# Plot the error
ggplot(
  data.frame(
    error = distance_classical_divide
  ), 
  aes(error)
) +
  geom_density()

summary(distance_classical_divide)