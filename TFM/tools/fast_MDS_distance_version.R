x = matrix(1:10^4, nrow = 100)
n = nrow(x)
p = 2
l = 20
s = 2
k = 2

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
  row.names(x) = as.character(1:nrow(x))
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
    list_matrix[[i_group]] = x[ind, ind]
  }
  

  able_to_do_mds = n/p <= l | p == 1
  
  
  # We can do MDS
  if(able_to_do_mds == TRUE){
    for (i_group in 1:p) {
      
      matrix_filter = list_matrix[[i_group]]
      distance_matrix = matrix_filter
      # MDS for each submatrix
      
      
      list_mds[[i_group]] = stats::cmdscale(
        d = distance_matrix, 
        k = s
      )
      
      
      # Subsample
      sample_size = sub_sample_size
      if(sample_size > length( row.names(matrix_filter ) ) ){
        sample_size = length( row.names(matrix_filter ) )
      }
      
      
      list_index[[i_group]] = sample(
        x = row.names(matrix_filter), 
        size = sample_size, 
        replace = FALSE
      )
      
      
      # Building indexes for M_align
      ind_M = which(row.names(x) %in% list_index[[i_group]])
      if(i_group == 1){
        indexes_M_align = ind_M
      }else{
        indexes_M_align = c(indexes_M_align, ind_M)
      }
    }
    
    distance_matrix_M = x[indexes_M_align, indexes_M_align]
    
    # M_align: MDS over x_M_align
    M_align = stats::cmdscale(
      d = distance_matrix_M, 
      k = s
    )
    
    
    # Global alignment
    for(i_group in 1:p){
      row_names = list_index[[i_group]]
      
      ind_M = which(row.names(M_align) %in% row_names)
      M_align_filter =  M_align[ind_M, ]
      
      di = list_mds[[i_group]]
      ind_mds = which(row.names( di ) %in% row_names)
      di_filter = di[ind_mds, ]
      
      # Alignment
      procrustes_result =  MCMCpack::procrustes(
        X = di_filter, #The matrix to be transformed
        Xstar = M_align_filter, # target matrix
        translation = TRUE, 
        dilation = TRUE
      )
      
      rotation_matrix = procrustes_result$R
      dilation = procrustes_result$s
      translation = procrustes_result$tt
      ones_vector = rep(1, nrow(di)) 
      translation_matrix = ones_vector %*% t(translation)
      
      
      tranformation_di = dilation * di %*% rotation_matrix + translation_matrix
      
      
      # Append
      if(i_group == 1){
        Z = tranformation_di
      } else{
        Z = rbind(
          Z,
          tranformation_di
        )
        
      }
    }
    
    row.names(Z) = row.names(x)
    
  }else{
    message("Recursive!!!")
    list_zi <- list()
    list_index <- list()
    
    for(i_group in 1:p){
      # Apply the algorithm
      list_zi[[i_group]] = fast_mds(
        x = list_matrix[[i_group]],
        n = nrow(list_matrix[[i_group]]),
        l = l,
        s = s,
        k = k,
        metric = metric 
      )
      
      #Take a subsample
      list_index[[i_group]] = sample(
        x = row.names( list_zi[[i_group]] ), 
        size = k * s, 
        replace = FALSE
      )
      message(paste0("End of getting z_", i_group))
      message(paste0("        Subsample is ", paste0( list_index[[i_group]], collapse = ',')))
      
      ind = which( row.names( list_zi[[i_group]] ) %in% list_index[[i_group]])
      
      
      message(paste0("        while bulding submatrix, its rownames are "), paste0(row.names(submatrix), collapse = ","))
      
      if(i_group == 1){
        indexes_M_align = ind 
      } else{
        indexes_M_align = rbind(
          indexes_M_align,
          ind
        )
      }
      message(paste0("        After iteration ", i_group, " x_M_align has ", nrow(x_M_align), " rows"))
    }
    message(paste0("        At the end x_M_align has ", nrow(x_M_align), " rows"))
    message(paste0("        At the end x_M_align has the following row names", paste0(row.names(x_M_align), collapse = ',' ) ))
    
    
    # M_align: MDS over x_M_align
    distance_matrix_M = x[indexes_M_align, indexes_M_align]
    
    
    M_align = stats::cmdscale(
      d = distance_matrix_M, 
      k = s
    )
    
    message(paste0("M_align for iterative has ", nrow(M_align)," rows"))
    message(paste0("x_M_align rows are ", paste0(row.names(x_M_align), collapse = ',' )))
    # Global alignment
    for(i_group in 1:p){
      row_names = list_index[[i_group]]
      message(paste0("selecting the following rows: ", paste0(row_names, collapse = ',')))
      ind_M = which(row.names(x_M_align) %in% row_names)
      M_align_filter =  M_align[ind_M, ]
      
      di = list_zi[[i_group]]
      ind_mds = which(row.names( di ) %in% row_names)
      di_filter = di[ind_mds, ]
      
      # Alignment
      message(paste0("       di_filter has: ", nrow(di_filter)))
      message(paste0("       M_align_filter has: ", nrow(M_align_filter)))
      procrustes_result =  MCMCpack::procrustes(
        X = di_filter, #The matrix to be transformed
        Xstar = M_align_filter, # target matrix
        translation = TRUE, 
        dilation = TRUE
      )
      
      rotation_matrix = procrustes_result$R
      dilation = procrustes_result$s
      translation = procrustes_result$tt
      ones_vector = rep(1, nrow(di)) 
      translation_matrix = ones_vector %*% t(translation)
      
      
      tranformation_di = dilation * di %*% rotation_matrix + translation_matrix
      
      
      # Append
      if(i_group == 1){
        Z = tranformation_di
      } else{
        Z = rbind(
          Z,
          tranformation_di
        )
        
      }
    }
    
    row.names(Z) = row.names(x)
    
  }
  
  return(Z)
}