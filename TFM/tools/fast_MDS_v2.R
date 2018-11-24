source("tools/load_libraries.R")

# Global variables
list_mds <<- list()
n_recursive_calls <<- 1
list_sample_points <<- list()
M_align <<- data.frame()


# This function decides whether it is possible to run classical mds or not
is_possible_to_run_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout = 60
){
  plan(multisession, gc = TRUE)
  is_possible_mds = TRUE
  
  # Calculate distance
  # message("Start calling the future")
  res_distance %<-% {
    daisy(
      x = x,
      metric = metric
    )
  }
  
  # message("Finish creation of future")
  is_able_to_resolve = NULL
  is_able_to_resolve <- withTimeout(
    {
      resolved(res_distance)
    },
    timeout = timeout,
    onTimeout = "silent"
  )
  
  is_mds_computed = is.null(is_able_to_resolve) == FALSE && 
    is_able_to_resolve == TRUE
 
  gc()
  return(is_mds_computed)
  
}


recursive_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout = 60,
  amplification_factor = 3
){
  
  # Amplify the number of coordinates by a factor
  number_coordinates_amplified = amplification_factor * number_coordinates
  
  # Check if it is possible to run MDS
  is_possible_mds = NULL
  is_possible_mds = is_possible_to_run_mds(
    x = x,
    number_coordinates = number_coordinates_amplified,
    metric = metric,
    timeout = timeout
  )
  
  # It is possible to run the MDS, the do it
  if(is_possible_mds == TRUE){
    message(paste0("Possible to run MDS with ", nrow(x), " rows"))
    message(paste0("Value of n_recursive_calls: ", n_recursive_calls))
    
    
    # Compute classical MDS
    
    # Calculate distance
    distance_matrix = daisy(
      x = x,
      metric = metric
    )
    
    # Applying MDS to the distance matrix
    classical_mds = cmdscale(
      d = distance_matrix, 
      eig = TRUE, 
      k = number_coordinates
    ) 
    
    output_mds = classical_mds$points
    
    # Save the output
    list_mds[[n_recursive_calls]] <<- output_mds
    
    # Take s points from the matrix randomly
    sampling_points = 0.9*nrow(x)
    samled_points = sample(
      x = row.names(x), 
      size = sampling_points, 
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
        timeout =  timeout,
        amplification_factor = 3
      )
    }
  }
}

fast_mds <- function(
  x,
  number_coordinates,
  metric,
  timeout = 60,
  amplification_factor = 3
){
  
  #Calling the recursive process
  recursive_mds(
    x = x,
    number_coordinates = number_coordinates,
    metric = metric,
    timeout = timeout,
    amplification_factor = amplification_factor
  )
  message("Finished the recursive MDS")
  
  # Once it is done, we stitch the solutions so that all the points have the
  # same system of coordinates. 
  
  # Perform MDS on M_align
  distance_M_align = daisy(
    x = M_align,
    metric = metric
  )
  
  message("Finished distance for M_align")
  
  mds_M_align = cmdscale(
    d = distance_M_align, 
    eig = TRUE, 
    k = number_coordinates
  )
  
  message("Finished MDS for M_align")
  
  mds_M_align = mds_M_align$points
    
  # We stitch the solutions
  total_groups = length(list_mds)
  mds_aligned = data.frame()
  for(i_group in 1:total_groups){
    # Filter the s points from D_i and M_align
    rows_to_select = list_sample_points[[i_group]]
    d_mds_i = list_mds[[i_group]]
    d_mds_i_filter = d_mds_i[rows_to_select,] 
    M_align_filter = mds_M_align[rows_to_select, ]
    
    # Solving procruster problem
    procusstes_result = procrustes(
      M_align_filter, 
      d_mds_i_filter
    )
    rotation_di = d_mds_i %*% procusstes_result$Q
    mds_aligned = rbind(
      mds_aligned,
      rotation_di
    )
  }
  
  
  # Reordering
  reording_permutation = match(row.names(x), row.names(mds_aligned))
  mds_aligned = mds_aligned[reording_permutation, ]
  
  
  return(mds_aligned)
  
}


###########################################################################
# Checks for fast_mds function
n_obs = 3*10^3
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



head(results_fast_mds, 100)


distance_x = daisy(
  x = x,
  metric = "euclidean"
)


mds_x = cmdscale(
  d = distance_x, 
  eig = TRUE, 
  k = 2
)


mds_x = mds_x$points


rotation_matrix = pracma::procrustes(
  mds_x,
  as.matrix(results_fast_mds)
)

rotated_results_fast_mds = as.matrix(results_fast_mds) %*% rotation_matrix$Q

df_rotated_results_fast_mds = as.data.frame(rotated_results_fast_mds)
df_mds_x = as.data.frame(mds_x)

df_rotated_results_fast_mds$source = 'fastMDS'
df_mds_x$source = 'ClassicalMDS'

all_results = rbind(
  df_rotated_results_fast_mds[1:10, ],
  df_mds_x[1:10,]
)

all_results %>% 
  ggplot(aes(x = V1, y = V2, group = source, color = source)) +
  geom_point() 

# Checks for is_possible_to_run_mds function
df_n_rows = seq(
  from = 10^3,
  to = 10^5,
  by = 10^3
)


df = data.frame(
  n_iter = rep(NA, length(df_n_rows)),
  is_possible_mds = rep(NA, length(df_n_rows))
)

i = 1
cont_it = TRUE
while(i <= length(df_n_rows) ){
  i_row = df_n_rows[i]
  # message("-------------------------------")
  # message(paste0("Nrow: ", i_row, " number of iteration ", i, " out of ", length(df_n_rows)))
  x = data.frame(
    x1 = rnorm(i_row),
    x2 = rnorm(i_row),
    x3 = rnorm(i_row),
    x4 = rnorm(i_row),
    x5 = rnorm(i_row),
    x6 = rnorm(i_row)
  )
  
  is_possible = is_possible_to_run_mds(
    x,
    number_coordinates = 2,
    metric = "euclidean",
    timeout = 5
  )
  
  df$n_iter[i] = i_row
  df$is_possible_mds[i] = is_possible
  
  cont_it = is_possible
  message(paste0("the result for this iteration is ", cont_it))
  message("----------------------------------------------------")
  i = i + 1
}




