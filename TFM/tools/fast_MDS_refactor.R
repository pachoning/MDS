fast_mds <- function(
  x,
  n,
  l,
  s,
  k,
  metric
){
  # Parametres inicials
  list_matrix = list()
  list_index = list()
  list_mds = list()
  list_mds_align = list()
  
  
  sub_sample_size = k * s
  
  
  # Division into p matrices
  # Puede ser que al hacer la particion, haya tantas matrices que k*s< nrow(x_i)
  # En este caso, volvemos a hacer un sampling
  p = ceiling(l/sub_sample_size)
  observations_division = sample(x = p, size = nrow(x), replace = TRUE)
  observations_division = sort(observations_division)
  min_sample_size = min(table(observations_division))
  
  while( min_sample_size < sub_sample_size && p > 1){
    p = p - 1 
    observations_division = sample(x = p, size = nrow(x), replace = TRUE)
    observations_division = sort(observations_division)
    min_sample_size = min(table(observations_division))
  }
  
  
  
  # Partition into p submatrices
  for(i_group in 1:p){
    ind = which(observations_division == i_group)
    list_matrix[[i_group]] = x[ind, ]
  }
  
  able_to_do_mds = n/p <= l | p == 1
  
  # We can do MDS
  if(able_to_do_mds == TRUE){
    for (i_group in 1:p) {
      
      matrix_filter = list_matrix[[i_group]]
      
      # MDS for each submatrix
      distance_matrix = daisy(
        x = matrix_filter,
        metric = metric
      )
      
      list_mds[[i_group]] = stats::cmdscale(
        d = distance_matrix, 
        k = s
      )
      
      #Take a subsample
      list_index[[i_group]] = sample(
        x = row.names( list_zi[[i_group]] ), 
        size = k * s, 
        replace = FALSE
      )
    }
    
  }else{
    message("Recursive!!!")
    list_zi <- list()
    list_index <- list()
    
    for(i_group in 1:p){
      # Apply the algorithm
      list_mds[[i_group]] = fast_mds(
        x = list_matrix[[i_group]],
        n = nrow(list_matrix[[i_group]]),
        l = l,
        s = s,
        k = k,
        metric = metric 
      )
      
      #Take a subsample
      list_index[[i_group]] = sample(
        x = row.names( list_matrix[[i_group]] ), 
        size = k * s, 
        replace = FALSE
      )
    }
    
    # Build x_M_align
    for(i_group in 1:p){
      
    }
    
  }
}