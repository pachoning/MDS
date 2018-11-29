source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
library(Ecdat)


recursive_divide_conquer<- function(
  x,
  n,
  l,
  p,
  s
){  
  # Initial parameters
  list_index <- list()
  
  # Partition of x in submatrices
  observations_division = sample(x = p, size = nrow(x), replace = TRUE)
  observations_division = sort(observations_division)
  min_sample_size = min(table(observations_division))
  
  while( min_sample_size < s && p > 1){
    p = p - 1 
    observations_division = sample(x = p, size = nrow(x), replace = TRUE)
    min_sample_size = min(table(observations_division))
  }
  
  
  
  # Partition into p submatrices
  for(i_group in 1:p){
    ind = which(observations_division == i_group)
    list_index[[i_group]] = row.names(x)[ind]
  }
  
  is_possible_to_run_divide = 2*n/p < l
  
  
  if(is_possible_to_run_divide == TRUE){
    message(paste0("Entering in the non-recursive part"))
    message(paste0("        n", n))
    message(paste0("        l", l))
    message(paste0("        2*n/p", 2*n/p))
    for(i_group in 1:p){
      
      #Selecting what to filter
      if(i_group == 1){
        index_to_filter = list_index[[i_group]]
      }else{
        index_to_filter = c(
          list_index[[i_group-1]],
          list_index[[i_group]]
        )
      }
      
      ind = which(row.names(x) %in% index_to_filter)
      matrix_filter = x[ind,]
      
      # MDS for each submatrix
      distance_matrix = daisy(
        x = matrix_filter,
        metric = metric
      )
      
      current_mds = stats::cmdscale(
        d = distance_matrix, 
        k = s
      )
      
      if(i_group == 1){
        cum_mds = current_mds
      }else{
        # Selecting the points from i-1 group in current_mds
        row_names_previous_group = list_index[[i_group-1]]
        ind_current_mds_previous_group = which(row.names(current_mds) %in% row_names_previous_group)
        current_mds_previous_group = current_mds[ind_current_mds_previous_group, ]
        
        # Selecting the points from i-1 group in cum_mds
        ind_cum_mds_previous_group = which(row.names(cum_mds) %in% row_names_previous_group)
        cum_mds_previous_group = cum_mds[ind_cum_mds_previous_group, ]
        
        # Selecting the points from i group in current_mds
        row_names_current_group = list_index[[i_group]]
        ind_current_mds_current_group = which(row.names(current_mds) %in% row_names_current_group)
        current_mds_current_group = current_mds[ind_current_mds_current_group, ]
        
        # Align the solutions
        procrustes_result =  MCMCpack::procrustes(
          X = current_mds_previous_group, #The matrix to be transformed
          Xstar = cum_mds_previous_group, # target matrix
          translation = TRUE, 
          dilation = TRUE
        )
        
        # Applying to the current group
        rotation_matrix = procrustes_result$R
        dilation = procrustes_result$s
        translation = procrustes_result$tt
        ones_vector = rep(1, nrow(current_mds_current_group)) 
        translation_matrix = ones_vector %*% t(translation)
        
        
        rotated_mds = dilation * current_mds_current_group %*% rotation_matrix + translation_matrix
        
        
        # Appending
        cum_mds = rbind(
          cum_mds,
          rotated_mds
        )
      }
      
      
    }
  } else{
    for(i_group in 1:p){
      # MDS of i
      ind = which(row.names(x) %in% list_index[[i_group]])
      message("entering in the recursive part")
      message(paste0("        n: ", nrow(x[ind, ])))
      message(paste0("        l: ", l))
      message(paste0("        p: ", p))
      message(paste0("        s: ", s))
      message(paste0("        2*n/p: ", 2*nrow(x[ind, ])/p))
      
      z = recursive_divide_conquer(
        x = x[ind, ],
        n = nrow(x[ind, ]),
        l = l,
        p = p,
        s = s
      )
      
      # MDS of i-1 and i
      if(i_group > 1){
        ind = which(row.names(x) %in% c( list_index[[i_group]], list_index[[i_group-1]]) )
        current_mds = recursive_divide_conquer(
          x = x[ind, ],
          n = nrow(x[ind, ]),
          l = l,
          p = p,
          s = s
        )
      }
      
      if(i_group == 1){
        cum_mds = z
      } else{
        # Selecting the points from i-1 group in cum_mds
        row_names_previous_group = list_index[[i_group - 1]]
        ind_cum_mds_previous_group = which(row.names(cum_mds) %in% row_names_previous_group)
        cum_mds_previous_group = cum_mds[ind_cum_mds_previous_group, ]
        
        # Selecting the points from i-1 group in the MDS that contains i-1 and i
        ind_current_mds_previous_group = which(row.names(current_mds) %in% row_names_previous_group)
        current_mds_previous_group = cum_mds[ind_current_mds_previous_group, ]
        
        
        # Selecting the points from i group in current_mds
        row_names_current_group = list_index[[i_group]]
        ind_current_mds_current_group = which(row.names(current_mds) %in% row_names_current_group)
        current_mds_current_group = current_mds[ind_current_mds_current_group, ]
        
        procrustes_result =  MCMCpack::procrustes(
          X = current_mds_previous_group, #The matrix to be transformed
          Xstar = cum_mds_previous_group, # target matrix
          translation = TRUE, 
          dilation = TRUE
        )
        
        # Applying to the current group
        rotation_matrix = procrustes_result$R
        dilation = procrustes_result$s
        translation = procrustes_result$tt
        ones_vector = rep(1, nrow(current_mds_current_group)) 
        translation_matrix = ones_vector %*% t(translation)
        
        rotated_mds = dilation * current_mds_current_group %*% rotation_matrix + translation_matrix
        
        
        # Appending
        cum_mds = rbind(
          cum_mds,
          rotated_mds
        )
      } 
    }
  }
  
  return(cum_mds)
}


x = BudgetFood %>% slice(1:300) %>% select(-sex, -town) 
n = nrow(x)
l = 180
p = 2
s = 2
metric = 'euclidean'

if(FALSE){
  results_classical_mds_budget = classical_mds(
    x = x,
    number_coordinates = 2,
    metric = metric
  )
}

2*n/p <= l

recursive_divide_conquer_result = recursive_divide_conquer(
    x = x,
    n = n,
    l = l,
    p = p,
    s = s
)

results_compare_divide_conquer_budget = compare_methods(
  mds_new_approach = recursive_divide_conquer_result,
  mds_classical = results_classical_mds_budget
)


head(recursive_divide_conquer_result, 8)
head(results_compare_divide_conquer_budget$mds_classical_transformed, 8)



