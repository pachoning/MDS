source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
library(Ecdat)

fast_roger <- function(
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
    message(paste0("Able to do MDS with p: ", p))
    message(paste0("       x has : ", nrow(x), " rows"))
    message(paste0("       Sample size of matrices: ", paste0(map_int(list_matrix, nrow), collapse = ',' ) ) )
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
  
      
      # Building x_M_align
      ind_M = which(row.names(x) %in% list_index[[i_group]])
      if(i_group == 1){
        x_M_align = x[ind_M, ]
      }else{
        x_M_align = rbind(
          x_M_align,
          x[ind_M, ]
        )
      }
      
    }
    
    # M_align: MDS over x_M_align
    distance_matrix_M = distance_matrix = daisy(
      x = x_M_align,
      metric = metric
    )
    
    
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

  }else{
    
    message("Recursive!!!")
    list_recursive <- list()
    for(i_group in 1:p){
      list_recursive[[i_group]] = fast_roger(
        x = list_matrix[[i_group]],
        n = nrow(list_matrix[[i_group]]),
        l = l,
        s = s,
        k = k,
        metric = metric 
      )
      
      if(i_group == 1){
        Z = list_recursive[[i_group]]
      }else{
        Z = rbind(
          Z,
          list_recursive[[i_group]]
        )
      }
    }
  }
  
  return(Z)
}



x = BudgetFood %>% slice(1:3000) %>% select(-sex, -town) 
n = nrow(x)
l = 150
s = 2
k = 3
metric = "euclidean"

n*k*s <= l^2


fast_r <- fast_roger(
  x = x,
  n = n,
  l = l,
  s = s,
  k = k,
  metric = "euclidean"
)

nrow(fast_r)

if(FALSE){
  results_classical_mds_budget = classical_mds(
    x = x,
    number_coordinates = 2,
    metric = metric
  )
}


results_compare_divide_conquer_budget = compare_methods(
  mds_new_approach = fast_r,
  mds_classical = results_classical_mds_budget
)


head(fast_r, 8)
head(results_compare_divide_conquer_budget$mds_classical_transformed, 8)

