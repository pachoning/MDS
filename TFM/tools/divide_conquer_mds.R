divide_conquer.mds <- function(
  x,
  l,
  s,
  metric
){
  
  # List positions
  ls_positions = list()
  list_eigenvalues = list()
  i_eigen = 1
  
  # Initial parameters
  p = ceiling(2*nrow(x)/l)
  groups = sample(x = p, size = nrow(x), replace = TRUE)
  groups = sort(groups)
  unique_group = unique(groups)
  total_groups = length(unique_group)
  
  for(k in 1:total_groups){
    # Getting the group that is being processed
    current_group = unique_group[k]
    x_positions_current_group = which(groups == current_group)
    ls_positions[[k]] = x_positions_current_group
    
    # Take the data in the following way:
    #   If it is the first iteration, take the data from first group
    #   else, take the data from k-1 and k groups
    if(k == 1){
      filter_rows_by_position = x_positions_current_group
      rows_processed = x_positions_current_group
    }else{
      
      rows_processed = c(
        rows_processed,
        x_positions_current_group
      )
      previous_group = unique_group[k-1]
      x_positions_previous_group =  which(groups == previous_group)
      
      # Rows to be filtered
      filter_rows_by_position = c(
        x_positions_previous_group,
        x_positions_current_group
      )
      
    }
    
    # Matrix to apply MDS
    submatrix_data = x[filter_rows_by_position, ]
    
    # Calculate distance
    distance_matrix = cluster::daisy(
      x = submatrix_data,
      metric = metric
    )
    
    # Applying MDS to the submatrix of data
    cmd_eig = stats::cmdscale(
      d = distance_matrix, 
      k = s,
      eig = TRUE
    )
    
    mds_iteration = cmd_eig$points
    if(p%%2 == 0){
      if(k%%2 == 0){
        list_eigenvalues[[i_eigen]] = cmd_eig$eig/nrow(submatrix_data)
        i_eigen = i_eigen + 1
      }
    }else{
      if(k %% 2 == 1){
        list_eigenvalues[[i_eigen]] = cmd_eig$eig/nrow(submatrix_data)
        i_eigen = i_eigen + 1
      }
    }

        row.names(mds_iteration) = row.names(submatrix_data)
    
    if(k == 1){
      # Define cum-MDS as MDS(1)
      cum_mds = mds_iteration
      
    }else{
      # Take the result of MDS(k-1) obtained with k-1 and k
      positions_previous_group_current_mds = row.names(mds_iteration) %in% row.names(x)[x_positions_previous_group]  
      positions_current_group_current_mds = row.names(mds_iteration) %in% row.names(x)[x_positions_current_group] 
      mds_previous = mds_iteration[positions_previous_group_current_mds,]  
      mds_current = mds_iteration[positions_current_group_current_mds,]
      
      # From cum-MDS take the result of group k-1
      positions_cum_sum_previous = which(row.names(cum_mds) %in% row.names(x)[x_positions_previous_group])
      cum_mds_previous = cum_mds[positions_cum_sum_previous, ] 
      
      
      # Apply Procrustes transformation
      procrustes_result =  MCMCpack::procrustes(
        X = mds_previous, #The matrix to be transformed
        Xstar = cum_mds_previous, # target matrix
        translation = TRUE, 
        dilation = TRUE
      )
      
      rotation_matrix = procrustes_result$R
      dilation = procrustes_result$s
      translation = procrustes_result$tt
      ones_vector = rep(1, nrow(mds_current)) 
      translation_matrix = ones_vector %*% t(translation)
      
      # Transform the data for the k-th group  
      cum_mds_current = dilation * mds_current %*% rotation_matrix + translation_matrix
      
      cum_mds = rbind(
        cum_mds,
        cum_mds_current
      )
      
    }
    
  }
  
  # Reordering
  reording_permutation = match(1:nrow(x), rows_processed)
  cum_mds = cum_mds[reording_permutation, ]
  
  
  return(
    list(
      points = cum_mds,
      eig = list_eigenvalues 
    )
  )
}