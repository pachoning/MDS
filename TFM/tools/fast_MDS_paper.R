source("tools/load_libraries.R")
n_times <<-  1
fast_mds_2 <- function(
  x,
  n,
  l,
  s,
  k,
  metric
){
  list_matrix = list()
  list_index = list()
  list_mds = list()
  
  sub_sample_size = k * s
  # Divide x in p submatrices
  # p = ceiling(l/sub_sample_size)
  p = ceiling(n/l)
  observations_division = sample(x = p, size = nrow(x), replace = TRUE)
  
  
  message(paste0("Calling algorithm with n: ", nrow(x)))
  message(paste0("Calculated p: ", p))
  message(paste0("n/p: ", n/p))
  
  
  # Partition into p matrices
  for(i_group in 1:p){
    ind = which(observations_division == i_group)
    list_index[[i_group]] = sample(x = row.names(x)[ind], size = sub_sample_size)
    list_matrix[[i_group]] = x[which(observations_division == i_group), ]
  }
  
  # If n/p<l, apply classical MDS
  if(n/p<l){
    for (i_group in 1:p) {
      distance_matrix = daisy(
        x = list_matrix[[i_group]],
        metric = metric
      )
      
      list_mds[[i_group]] = stats::cmdscale(
        d = distance_matrix, 
        k = s
      )
    }
  }else{
    message(paste0("calling recursively with n: ", n))
  # Call Fast MDS
    for(i_group in 1:p){
      fast_mds_2(
        x = list_matrix[[i_group]],
        n = nrow(list_matrix[[i_group]]),
        l = l,
        s = s,
        k = k,
        metric = metric 
      )
    }
  }
  message(paste0("after if else with n: ", nrow(x)))
  
  # Selecting s points from p submatrices
  for(i_group in 1:p){
    # Sampling s * k points from x
    ind = row.names(x) %in% list_index[[i_group]]
    sub_matrix =  x[ind, ]
    if(i_group == 1){
      x_M_align =sub_matrix
    }else{
      x_M_align = rbind(
        x_M_align,
        sub_matrix
      )
    }
  }
  
  # Getting M_align
  # Calculate distance
  distance_matrix = daisy(
    x = x_M_align,
    metric = metric
  )
  
  # M_alignm
  M_align =  stats::cmdscale(
    d = distance_matrix, 
    k = s
  ) 
  
  # Alignement
  for(i_group in 1:p){
    di = list_mds[[i_group]]
    ind_di = row.names(di) %in% list_index[[i_group]]
    di_filter = di[ind_di,]
    
    ind_Mi = row.names(M_align) %in% list_index[[i_group]]
    Mi_filter = M_align[ind_Mi, ]
    procrustes_result =  MCMCpack::procrustes(
      X = di_filter, #The matrix to be transformed
      Xstar = Mi_filter, # target matrix
      translation = TRUE, 
      dilation = TRUE
    )
    
    rotation_matrix = procrustes_result$R
    dilation = procrustes_result$s
    translation = procrustes_result$tt
    ones_vector = rep(1, nrow(di)) 
    translation_matrix = ones_vector %*% t(translation)
    
    rotation_matrix = procrustes_result$R
    dilation = procrustes_result$s
    translation = procrustes_result$tt
    ones_vector = rep(1, nrow(di)) 
    translation_matrix = ones_vector %*% t(translation)
    
    
    # Transforming the data for the k-th group  
    tranformation_di = dilation * di %*% rotation_matrix + translation_matrix
  
    if(i_group == 1){
      Z = tranformation_di
    }else{
      Z = rbind(
        Z,
        tranformation_di
      )
    }
  }
  
  message(paste0("rows of Z:", nrow(Z)))
  permutation_to_order = match(row.names(x), row.names(Z))
  return(Z[permutation_to_order, ])
}


n = 2000
l = 1000
s = 2
k = 3
metric = "euclidean"
# Toy matrix
x = matrix(rnorm(5*n), nrow = n)
row.names(x) = 1:n


fast_mds_sol <- fast_mds_2(
  x = x,
  n = nrow(x),
  l = l,
  s = 2,
  k = 3,
  metric = "euclidean"
)

