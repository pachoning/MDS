# This function decides whether it is possible to compute distance matrix
is_possible_to_calculate_distance_matrix <- function(
  x,
  metric,
  timeout
  
){
  
  
  
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
  
  
  if(is_distance_computed == FALSE){
    distance_matrix = NULL
  }else{
    distance_matrix = res_distance
  }
  
  list_to_return = list(
    is_distance_computed = is_distance_computed,
    distance_matrix = distance_matrix
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
  plan(multisession, gc = TRUE)
  
  # Check if we can compute distance matrix
  result_distance_computation = is_possible_to_calculate_distance_matrix(
    x = x,
    metric = metric,
    timeout = timeout
  )
  

  able_to_compute_distance = result_distance_computation$is_distance_computed
  able_to_compute_components_mds = able_to_compute_distance
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
  
  # Store the rows of the main data frame
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
    sampled_points = sample(
      x = row.names(x), 
      size = n_sampling_points, 
      replace = FALSE
    )
    
    list_sample_points[[n_recursive_calls]] <<- sampled_points
    
    # Append the points no M_align
    ind_position = row.names(x) %in% sampled_points
    if(is_M_align_empty == TRUE){
      M_align <<- x[ind_position, ]
      row.names(M_align) = row.names(x)[ind_position]
      is_M_align_empty <<- FALSE
    }else{
      M_to_append = x[ind_position, ]
      row.names(M_to_append) = row.names(x)[ind_position]
      M_align <<- rbind(
        M_align,
        M_to_append
      )
      
    }
    

    
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
  timeout

){
  
  # Global variables
  list_mds <<- list()
  n_recursive_calls <<- 1
  list_sample_points <<- list()
  is_M_align_empty <<- TRUE
  M_align <<- NA
  n_obs_main_matrix <<- NA
  n_sampling_points <<- NA
  
  
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
    rows_names_to_select = list_sample_points[[i_group]]
    d_mds_i = list_mds[[i_group]]
    rows_to_select_d = which(row.names(d_mds_i) %in% rows_names_to_select)
    rows_to_select_M = which(row.names(mds_M_align) %in% rows_names_to_select)
    d_mds_i_filter = d_mds_i[rows_to_select_d,] 
    M_align_filter = mds_M_align[rows_to_select_M, ]
    
    # Solving procruster problem
    procrustes_result =  MCMCpack::procrustes(
      X = d_mds_i_filter, #The matrix to be transformed
      Xstar = M_align_filter, # target matrix
      translation = TRUE, 
      dilation = TRUE
    )
    
    rotation_matrix = procrustes_result$R
    dilation = procrustes_result$s
    translation = procrustes_result$tt
    ones_vector = rep(1, nrow(d_mds_i)) 
    translation_matrix = ones_vector %*% t(translation)
    
    # Transforming the data for the k-th group  
    tranformation_di = dilation * d_mds_i %*% rotation_matrix + translation_matrix
    
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