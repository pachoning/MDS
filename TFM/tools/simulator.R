get.mean.eigenvalues <- function(
  list_eigenvalues,
  n_eigenvalues
){

  if( length(list_eigenvalues) <= 1 & 
      (
        is.null(list_eigenvalues) == TRUE || is.na(list_eigenvalues) == TRUE
      )
  ){
    vector_mean_eigen = NA
  }else{
    if( is.list(list_eigenvalues) == FALSE ){
      list_eigenvalues = list(list_eigenvalues)
    }
    
    vector_eigen = rapply(
      list_eigenvalues, 
      function(x, th = n_eigenvalues) x[1:th], 
      how = "unlist"
    )
    
    length_vector_eigen = length(vector_eigen)
    vector_mean_eigen = rep(NA, n_eigenvalues)
    vector_index = 1:length_vector_eigen
    for(i in 0:(n_eigenvalues-1)){
      ind = which(vector_index%%n_eigenvalues == i) 
      if(i == 0){
        j = n_eigenvalues
      }else{
        j = i
      }
      vector_mean_eigen[j] = mean(vector_eigen[ind], na.rm = TRUE)
    }
    
  }

  return(vector_mean_eigen)
}


n.dimensions <- function(
  list_eigenvalues,
  threshold_main_dimensions
){
  if( length(list_eigenvalues) <= 1 & 
      (
        is.null(list_eigenvalues) == TRUE || is.na(list_eigenvalues) == TRUE
      )
  ){
    above_threshold = NA
  }else{
    if( is.list(list_eigenvalues) == FALSE ){
      list_eigenvalues = list(list_eigenvalues)
    }
    
    above_threshold = rapply(
      list_eigenvalues, 
      function(x, th = threshold_main_dimensions) min(which(cumsum(x)/sum(x) > th)), 
      how = "unlist"
    )
    
    above_threshold = max(above_threshold)
  }
  
  return(above_threshold)
}


determine.depth <- function(this,thisdepth=0){
  if( !is.list(this) ){
    return(thisdepth)
  }else{
    return(
      max(
        unlist(
          lapply(
            this,
            determine.depth,
            thisdepth=thisdepth+1
          )
        )
      )
    )    
  }
}


# Function to generate the data
build.data <- function(
  sample_size,
  data_dimension,
  main_dimensions_vector
){
  
  # Build lambda vecotr (dimesionality)
  lambda_vector = rep(1, data_dimension)
  
  real_data_dimension = length(main_dimensions_vector)
  if( real_data_dimension > 0){
    if(real_data_dimension > data_dimension){
      stop("main_dimensions_vector is greater than data_dimension")
    }
    lambda_vector[1:real_data_dimension] = main_dimensions_vector
  } 
  
  # Buil data
  x = matrix(
    rnorm(
      n = sample_size*data_dimension
    ),
    ncol = data_dimension,
    nrow = sample_size
  )
  
  row.names(x) = as.character(1:sample_size)
  
  x = x %*% diag(lambda_vector)
  return(x)
}


aggregator.mds <- function(
  x,
  l,
  s, 
  k,
  metric,
  method_wanted
){
  
  starting_time = proc.time()
  if( method_wanted == 'divide_conquer' ){
    message("Performing divide and conquer mds")
    result_mds = divide_conquer.mds(
      x = x,
      l = l,
      s = s,
      metric = metric
    )
  }else if( method_wanted == 'fast' ){
    message("Performing fast mds")
    result_mds = fast.mds(
      x = x,
      n = nrow(x),
      l = l,
      s = s,
      k = k,
      metric = metric
    )
    
  }else if( method_wanted == 'classical' ){
    message("Performing classical mds")
    result_mds = classical.mds(
      x = x,
      s = s,
      metric = metric
    )
  }else if( method_wanted == 'gower'){
    message("Performing Gower mds")
    result_mds = gower.interpolation.mds(
      x = x,
      l = l,
      s = s,
      metric = metric 
    )
  }else{
    stop( "invalid value for method_wanted variable" )
  }
  
  diff_time = proc.time() - starting_time 
  elapsed_time_method = round(diff_time[3], 4)
  return(
    list(
      points = result_mds$points,
      eig = result_mds$eig,
      elapsed_time = elapsed_time_method
    )
  )
}

canonical.corr <- function(
  x,
  y
  
){
  can_corr <- CCA::cc(x, y)
  can_corr_result = min(can_corr$cor)
  return(can_corr_result)
}


corr.groups.procrustes <- function(
  x_to_be_transformed,
  x_target
){
  
  if(is.matrix(x_to_be_transformed) == FALSE) x_to_be_transformed = as.matrix(x_to_be_transformed)
  if(is.matrix(x_target) == FALSE) x_target = as.matrix(x_target)
  
  procrustes_analysis =  MCMCpack::procrustes(
    X = x_to_be_transformed, #The matrix to be transformed
    Xstar = x_target, # target matrix
    translation = TRUE, 
    dilation = TRUE
  )
  
  return(
    diag(
      cor(
        x_target, 
        procrustes_analysis$X.new
      )
    )
  )
}



# Generation of data
do.magic <- function(
  sample_size,
  data_dimension, # Number of columns
  main_dimensions_vector, # Number of real dimension
  n_primary_dimensions,
  l,
  k,
  metric,
  compute_divide_conquer_mds,
  compute_fast_mds,
  compute_gower_mds,
  compute_classical_mds,
  max_sample_size_classical,
  threshold_main_dimensions,
  n_eigenvalues,
  n_cols_procrustes_noise
){
  # Dimensions to do procrustes
  if(n_primary_dimensions == 0){
   n_dimensions_procrustes = n_cols_procrustes_noise 
  }else{
    n_dimensions_procrustes = n_primary_dimensions 
  }
  
  # Build the data 
  x <- build.data(
    sample_size = sample_size,
    data_dimension = data_dimension,
    main_dimensions_vector = main_dimensions_vector
  )


  # Run divide and conquer in case it is needed
  divide_conquer_points = NA
  divide_conquer_eig = NA
  divide_conquer_eig_subsample = NA
  divide_conquer_n_dimensions = NA
  divide_conquer_elapsed_time = NA
  divide_conquer_corr_matrix = NA
  
  
  if( compute_divide_conquer_mds == TRUE ){
    divide_conquer_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'divide_conquer'
    )
    
    divide_conquer_points = divide_conquer_mds$points
    divide_conquer_eig = divide_conquer_mds$eig
    divide_conquer_eig_subsample = get.mean.eigenvalues(
      list_eigenvalues = divide_conquer_eig,
      n_eigenvalues =n_eigenvalues
    )
    divide_conquer_n_dimensions = n.dimensions(
      list_eigenvalues = divide_conquer_eig,
      threshold_main_dimensions = threshold_main_dimensions
    )
    divide_conquer_elapsed_time = as.numeric(divide_conquer_mds$elapsed_time)
    divide_conquer_corr_matrix = corr.groups.procrustes(
      x_to_be_transformed = x[, 1:n_dimensions_procrustes],
      x_target = divide_conquer_points[, 1:n_dimensions_procrustes]
    )
    
  }

  # Run fast in case it is needed
  fast_points = NA
  fast_eig = NA
  fast_eig_subsample = NA
  fast_n_dimensions = NA
  fast_elapsed_time = NA
  fast_corr_matrix = NA
 
  
  if( compute_fast_mds == TRUE ){
    fast_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'fast'
    )
    
    fast_points = fast_mds$points
    fast_eig = fast_mds$eig
    fast_eig_subsample = get.mean.eigenvalues(
      list_eigenvalues = fast_eig,
      n_eigenvalues =n_eigenvalues
    )
    fast_n_dimensions = n.dimensions(
      list_eigenvalues = fast_eig,
      threshold_main_dimensions = threshold_main_dimensions
    )
    fast_elapsed_time = as.numeric(fast_mds$elapsed_time)
    fast_corr_matrix = corr.groups.procrustes(
      x_to_be_transformed = x[, 1:n_dimensions_procrustes],
      x_target = fast_points[, 1:n_dimensions_procrustes]
    )
  }
  
  # Run Gower
  gower_points = NA
  gower_eig = NA
  gower_eig_subsample = NA
  gower_n_dimensions = NA
  gower_elapsed_time = NA
  gower_corr_matrix = NA
  
  if( compute_gower_mds == TRUE ){
    gower_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'gower'
    )
    
    gower_points = gower_mds$points
    gower_eig = gower_mds$eig
    gower_eig_subsample = get.mean.eigenvalues(
      list_eigenvalues = gower_eig,
      n_eigenvalues = n_eigenvalues
    )
    gower_n_dimensions = n.dimensions(
      list_eigenvalues = gower_eig,
      threshold_main_dimensions = threshold_main_dimensions
    )
    gower_elapsed_time = as.numeric(gower_mds$elapsed_time)
    gower_corr_matrix = corr.groups.procrustes(
      x_to_be_transformed = x[, 1:n_dimensions_procrustes],
      x_target = gower_points[, 1:n_dimensions_procrustes]
    )
    
  }


  # Run the classical in case it is needed
  classical_points = NA
  classical_eig = NA
  classical_eig_subsample = NA
  classical_n_dimensions = NA
  classical_elapsed_time = NA
  classical_corr_matrix = NA
  
  if( compute_classical_mds == TRUE ){
    classical_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'classical'
    )
    
    classical_points = classical_mds$points
    classical_eig = classical_mds$eig
    classical_eig_subsample = get.mean.eigenvalues(
      list_eigenvalues = classical_eig,
      n_eigenvalues = n_eigenvalues
    )
    classical_n_dimensions = n.dimensions(
      list_eigenvalues = classical_eig,
      threshold_main_dimensions = threshold_main_dimensions
    )
    classical_elapsed_time = as.numeric(classical_mds$elapsed_time)
    classical_corr_matrix = corr.groups.procrustes(
      x_to_be_transformed = x[, 1:n_dimensions_procrustes],
      x_target = classical_points[, 1:n_dimensions_procrustes]
    )
    
  }

  list_to_return = list(
    # Inputs
    x = x,
    sample_size = sample_size,
    data_dimension = data_dimension,
    main_dimensions_vector = main_dimensions_vector,
    l = l,
    k = k,
    metric = metric,
    
    # Output for divide and conquer
    divide_conquer_points = divide_conquer_points,
    divide_conquer_eig_subsample = divide_conquer_eig_subsample,
    divide_conquer_n_dimensions = divide_conquer_n_dimensions,
    divide_conquer_elapsed_time = divide_conquer_elapsed_time,
    divide_conquer_corr_matrix = divide_conquer_corr_matrix,
    
    # Output for fast
    fast_points = fast_points,
    fast_eig_subsample = fast_eig_subsample,
    fast_n_dimensions = fast_n_dimensions,
    fast_elapsed_time = fast_elapsed_time,
    fast_corr_matrix = fast_corr_matrix,
    
    # Output form gower
    gower_points = gower_points,
    gower_eig_subsample = gower_eig_subsample,
    gower_n_dimensions = gower_n_dimensions,
    gower_elapsed_time = gower_elapsed_time,
    gower_corr_matrix = gower_corr_matrix,
    
    
    # Output for classical
    classical_points = classical_points,
    classical_eig_subsample = classical_eig_subsample,
    classical_n_dimensions = classical_n_dimensions,
    classical_elapsed_time = classical_elapsed_time,
    classical_corr_matrix = classical_corr_matrix
  )
  
  return(list_to_return)
}
