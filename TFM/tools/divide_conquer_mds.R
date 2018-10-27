divide_conquer_mds <- function(
  x,
  groups,
  number_coordinates,
  metric = "gower",
  ...
){
  # Loading libraries
  if( require(cluster) == FALSE ){
    install.packages('cluster')
    library(cluster)
  }
  
  if( require(pracma) == FALSE ){
    install.packages('pracma')
    library(pracma)
  }
  
  
  # Initial parameters
  unique_group = unique(groups)
  total_groups = length(unique_group)
  error = c()
  
  for(k in 1:total_groups){
    # Getting the group that is being processed
    current_group = unique_group[k]
    positions_current_group = which(groups == current_group)
    total_elements_current_group = length(positions_current_group)
    
    # Take the data in the following way:
    #   If it is the first iteration, take the data from fist group
    #   else, take the data from k-1 and k groups
    if(k == 1){
      filter_rows_by_position = positions_current_group
    }else{
      previous_group = unique_group[k-1]
      positions_previous_group =  which(groups == previous_group)
      total_elements_previous_group = length(positions_previous_group)
      
      # Registers to be filtered
      filter_rows_by_position = c(
        positions_previous_group,
        positions_current_group
      )
      
      total_elements_mds = length(filter_rows_by_position)
    }
    message(paste0("Iteration number ", k, " out of ", total_groups))
    
    # Matrix to apply MDS
    submatrix_data = x[filter_rows_by_position, ]
    
    # Calculate distance
    distance_matrix = daisy(
      submatrix_data,
      metric = metric
    )
    
    # Applying MDS to the submatrix of data
    mds_iteration = cmdscale(
      distance_matrix, 
      eig = TRUE, 
      k = number_coordinates
    ) 
    
    mds_iteration = mds_iteration$points
    
    
    if(k == 1){
      # Define cum-MDS as MDS(1)
      cum_mds = mds_iteration
      positions_cum_sum = as.vector(rep(current_group, total_elements_current_group))
      
    }else{
      # Take the result of MDS(k-1) obtained with k-1 and k
      positions_previous_mds_iteration = 1:total_elements_previous_group
      positins_current_mds_iteration = (total_elements_previous_group+1):total_elements_mds
      
      mds_previous = mds_iteration[positions_previous_mds_iteration,]  
      mds_current = mds_iteration[positins_current_mds_iteration,]
      
      # From cum-MDS take the result of group k-1
      positions_cum_sum_previous = which(positions_cum_sum == previous_group)
      cum_mds_previous = cum_mds[positions_cum_sum_previous, ] 
      
      
      # Applying Procrustes transformation
      procusstes_result = procrustes(
        cum_mds_previous, 
        mds_previous
      )
      
      error[k-1] = procusstes_result$d
      
      translation_matrix = procusstes_result$Q
      
      # Rotating the data for the k-th group  
      cum_mds_current = mds_current %*% procusstes_result$Q
      cum_mds = rbind(
        cum_mds,
        cum_mds_current
      )
      
      new_indexes = as.vector(rep(current_group, total_elements_current_group))
      
      positions_cum_sum = c(
        positions_cum_sum,
        new_indexes
      )
      
    }
    
  }
  # Reordering
  reording_permutation = order(row.names(cum_mds), row.names(x) )
  cum_mds = cum_mds[reording_permutation, ]
  
  return(
    list(
      mds = cum_mds,
      error = error
    )
  )
}