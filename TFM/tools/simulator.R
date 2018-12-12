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
  
  procrustes_analysis =  MCMCpack::procrustes(
    X = x_to_be_transformed, #The matrix to be transformed
    Xstar = x_target, # target matrix
    translation = TRUE, 
    dilation = TRUE
  )
  
  return(
    min(
      diag(
        cor(
          x_target, 
          procrustes_analysis$X.new
        )
      )
    )
  )
}







# Generation of data
do.magic <- function(
  sample_size,
  data_dimension, # Number of columns
  main_dimensions_vector, # Number of real dimension
  l,
  k,
  metric,
  compute_divide_conquer_mds,
  compute_fast_mds,
  compute_classical_mds,
  max_sample_size_classical
){

  # Build the data 
  x <- build.data(
    sample_size = sample_size,
    data_dimension = data_dimension,
    main_dimensions_vector = main_dimensions_vector
  )


  # Run divide and conquer in case it is needed
  divide_conquer_points = NA
  divide_conquer_eig = NA
  divide_conquer_elapsed_time = NA
  divide_conquer_canonical_corr = NA
  divide_conquer_correlation_data = NA

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
    divide_conquer_elapsed_time = as.numeric(divide_conquer_mds$elapsed_time)
    
    # Canonical correlation between data and divide and conquer
    divide_conquer_canonical_corr = canonical.corr(
      x = x,
      y = divide_conquer_points
    )
    
    # Correlation between data and divide and conquer
    divide_conquer_correlation_data = corr.groups.procrustes(
      x_to_be_transformed = x,
      x_target = divide_conquer_points
    )
    
  }

  # Run fast in case it is needed
  fast_points = NA
  fast_eig = NA
  fast_elapsed_time = NA
  fast_canonical_corr = NA
  fast_correlation_data = NA
  
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
    
    # Unlist the eigenvalues
    fast_eig = fast_mds$eig
    fast_elapsed_time = as.numeric(fast_mds$elapsed_time)
    
    # Canonical correlation between data and fast
    fast_canonical_corr = canonical.corr(
      x = x,
      y = fast_points
    )
    
    # Correlation between data and fast
    fast_correlation_data = corr.groups.procrustes(
      x_to_be_transformed = x,
      x_target = fast_points
    )
  }


  # Run the classical in case it is needed
  classical_points = NA
  classical_eig = NA
  classical_elapsed_time = NA
  classical_canonical_corr = NA
  classical_correlation_data = NA
  
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
    classical_elapsed_time = as.numeric(classical_mds$elapsed_time)
    
    # Canonical correlation between data and classical MDS
    classical_canonical_corr = canonical.corr(
      x = x,
      y = classical_points
    )
    
    # Correlation between data and classical MDS
    classical_correlation_data = corr.groups.procrustes(
      x_to_be_transformed = x,
      x_target = classical_points
    )
  }
  
  # Correlations between divide and conquer and MDS classical
  divide_conquer_corr_classical_mds = NA
  if( compute_classical_mds == TRUE && compute_divide_conquer_mds ==TRUE ){
    # Correlations between classical MDS and divide and conquer
    divide_conquer_corr_classical_mds = corr.groups.procrustes(
      x_to_be_transformed = classical_points,
      x_target = divide_conquer_points
    )
  }
  
  # Correlations between classical MDS and fast
  fast_corr_mds_classical = NA
  if( compute_classical_mds == TRUE && compute_fast_mds ==TRUE ){
    fast_corr_mds_classical = corr.groups.procrustes(
      x_to_be_transformed = classical_points,
      x_target = fast_points
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
    divide_conquer_eig = divide_conquer_eig,
    divide_conquer_elapsed_time = divide_conquer_elapsed_time,
    divide_conquer_canonical_corr = divide_conquer_canonical_corr,
    divide_conquer_correlation_data = divide_conquer_correlation_data,
    divide_conquer_corr_classical_mds = divide_conquer_corr_classical_mds,
    
    # Output for fast
    fast_points = fast_points,
    fast_eig = fast_eig,
    fast_elapsed_time = fast_elapsed_time,
    fast_canonical_corr = fast_canonical_corr,
    fast_correlation_data = fast_correlation_data,
    fast_corr_mds_classical = fast_corr_mds_classical,
    
    # Output for classical
    classical_points = classical_points,
    classical_eig = classical_eig,
    classical_elapsed_time = classical_elapsed_time,
    classical_canonical_corr = classical_canonical_corr,
    classical_correlation_data = classical_correlation_data
  )
  
  return(list_to_return)
}
