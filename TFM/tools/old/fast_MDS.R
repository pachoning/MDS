source("tools/load_libraries.R")

# Global variables
list_mds <<- list()
list_positions <<- list()
n_recursive_calls <<- 0


# This function decides whether it is possible to run classical mds or not
is_possible_to_run_mds <- function(
  x,
  number_coordinates,
  metric = "gower",
  timeout = 60
){
  
  is_possible_mds = TRUE
  
  # Calculate distance
  res_distance = NULL
  res_distance <- withTimeout(
    {
      distance_matrix = daisy(
        x,
        metric = metric
      )
    }, 
    timeout = timeout, 
    onTimeout = "silent"
  )
  
  # Not able to compute the distance matrix
  if(is.null(res_distance) == TRUE) is_possible_mds = FALSE
  
  
  # Computing MDS. It is calculated if it was able to compute the distance
  # matrix. Otherwise, this part is not executed
  # if(is_possible_mds == TRUE){
  #   res_mds = NULL
  #   
  #   res_mds <- withTimeout(
  #     {
  #       mds_iteration = cmdscale(
  #         distance_matrix, 
  #         eig = TRUE, 
  #         k = number_coordinates
  #       ) 
  #     }, 
  #     timeout = timeout, 
  #     onTimeout = "silent"
  #   )
  #   
  #   if(is.null(res_mds) == TRUE) is_possible_mds = FALSE
  # }
  
  return(is_possible_mds)
  
}


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
    timeout = 1
  )
  
  df$n_iter[i] = i_row
  df$is_possible_mds[i] = is_possible
  message(paste0("iteration ", i, " out of ", length(df_n_rows)))
  cont_it = is_possible
  message(paste0("the result for this iteration is ", cont_it))
  i = i + 1
}

head(df)

# Remember to multiple the number of dimensions by a factor of 3



?evalWithTimeout
?tryCatch
?evalWithTimeout
