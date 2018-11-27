source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
library(Ecdat)

is_computed_first_time = FALSE
Z_mds <<- NA


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
  p = ceiling(l/sub_sample_size)
  # p = ceiling(n/l)
  # p = 2
  observations_division = sample(x = p, size = nrow(x), replace = TRUE)
  
  message(paste0("Calling algorithm with n: ", nrow(x)))
  message(paste0("Calculated p: ", p))
  message(paste0("n/p: ", n/p))
  
  
  # Partition into p matrices
  for(i_group in 1:p){
    ind = which(observations_division == i_group)
    list_index[[i_group]] = sample(x = row.names(x)[ind], size = sub_sample_size, replace = FALSE)
    list_matrix[[i_group]] = x[ind, ]
  }
  
  # If n/p<l, apply classical MDS
  if(n/p<=l){
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
  
  
  if(length(list_mds) > 0){
    # Selecting points from p submatrices
    for(i_group in 1:p){
      # Sampling s * k points from x
      ind = row.names(x) %in% list_index[[i_group]]
      sub_matrix =  x[ind, ]
      if(i_group == 1){
        x_M_align = sub_matrix
      }else{
        x_M_align = rbind(
          x_M_align,
          sub_matrix
        )
      }
    }
    
    # Getting M_align
    # Calculate distance
    distance_matrix_M = daisy(
      x = x_M_align,
      metric = metric
    )
  
    # M_alignm
    M_align =  stats::cmdscale(
      d = distance_matrix_M, 
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
    Z = Z[permutation_to_order, ]
    if(is_computed_first_time == FALSE){
      is_computed_first_time <<- TRUE
      Z_mds <<- Z
      
    }else{
      Z_mds <<- rbind(
        Z_mds,
        Z
      )
    }
  }
}


x = BudgetFood %>% slice(1:3000) %>% select(-sex, -town) 
metric = "euclidean"
nrow(x)


fast_mds_sol <- fast_mds_2(
  x = x,
  n = nrow(x),
  l = 1000,
  s = 2,
  k = 3,
  metric = "euclidean"
)

nrow(Z_mds)
head(Z_mds)

if(FALSE){
  results_classical_mds_budget = classical_mds(
    x = x,
    number_coordinates = 2,
    metric = metric
  )
}


results_compare_divide_conquer_budget = compare_methods(
  mds_new_approach = Z_mds,
  mds_classical = results_classical_mds_budget
)


head(Z_mds, 8)
head(results_compare_divide_conquer_budget$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer_budget$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)




