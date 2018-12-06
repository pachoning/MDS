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
  method_wanted,
  sample_size_classical
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
    if(is.null(sample_size_classical) == TRUE){
      sample_size = nrow(x)
    }else{
      sample_size = sample_size_classical
    }
    
    rows_filter = 1:sample_size
    rows_filter = sort(rows_filter)
    
    x_filter = x[rows_filter, ]
    result_mds = classical.mds(
      x = x_filter,
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
  sample_size_classical
){

  # Build the data 
  x <- build.data(
    sample_size = sample_size,
    data_dimension = data_dimension,
    main_dimensions_vector = main_dimensions_vector
  )


  # Run divide and conquer in case it is needed
  divide_conquer_points = NULL
  divide_conquer_eig = NULL
  divide_conquer_elapsed_time = NULL


  if( compute_divide_conquer_mds == TRUE ){
    divide_conquer_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'divide_conquer',
      sample_size_classical = sample_size_classical
    )
    
    divide_conquer_points = divide_conquer_mds$points
    divide_conquer_eig = divide_conquer_mds$eig
    divide_conquer_elapsed_time = as.numeric(divide_conquer_mds$elapsed_time)
  }

  # Run fast in case it is needed
  fast_points = NULL
  fast_eig = NULL
  fast_elapsed_time = NULL
  
  if( compute_fast_mds == TRUE ){
    fast_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'fast',
      sample_size_classical = sample_size_classical
    )
    
    fast_points = fast_mds$points
    
    # Unlist the eigenvalues
    fast_eig = fast_mds$eig
    
    depth_list = determine.depth(fast_eig)
    
    if(depth_list > 1){
      for(i in 1:(depth_list-1)){
        fast_eig = unlist(fast_eig, recursive = FALSE)
      }
    }
    
    fast_elapsed_time = as.numeric(fast_mds$elapsed_time)
  }


  # Run the classical in case it is needed
  classical_points = NULL
  classical_eig = NULL
  classical_elapsed_time = NULL
  if( compute_classical_mds== TRUE ){
    classical_mds = aggregator.mds(
      x = x,
      l = l,
      s = data_dimension, 
      k = k,
      metric = metric,
      method_wanted = 'classical',
      sample_size_classical = sample_size_classical
    )
    
    classical_points = classical_mds$points
    classical_eig = classical_mds$eig
    classical_elapsed_time = as.numeric(classical_mds$elapsed_time)
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
    sample_size_classical = sample_size_classical,
    
    # Output for divide and conquer
    divide_conquer_points = divide_conquer_points,
    divide_conquer_eig = divide_conquer_eig,
    divide_conquer_elapsed_time = divide_conquer_elapsed_time,
    
    
    # Output for fast
    fast_points = fast_points,
    fast_eig = fast_eig,
    fast_elapsed_time = fast_elapsed_time,
    
    
    # Output for classical
    classical_points = classical_points,
    classical_eig = classical_eig,
    classical_elapsed_time = classical_elapsed_time
  )
  
  return(list_to_return)
}
