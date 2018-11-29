source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
library(Ecdat)

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
  p = ceiling(l/sub_sample_size)
  
  message("Beginnig of algorithm: ")
  message(paste0("        Calling algorithm with n: ", nrow(x)))
  message(paste0("        Calculated p: ", p))
  message(paste0("        n/p: ", n/p))
  
  # Divide x in p submatrices
  message(paste0("        value of size(1) is  ", nrow(x)))
  observations_division = sample(x = p, size = nrow(x), replace = TRUE)
  observations_division = sort(observations_division)
  
  
  # Partition into p matrices
  for(i_group in 1:p){
    ind = which(observations_division == i_group)
    ind = sort(ind)
    
    # Get x1,..., xp 
    list_matrix[[i_group]] = x[ind, ]
    
    # Take the rows of matrice i for i in 1:p
    n_sample = sub_sample_size
    if( n_sample > nrow( x[ind, ]) ) n_sample = nrow( row.names(x)[ind]) 

    message(paste0("        value of size(2) is  ", n_sample))
    
    list_index[[i_group]] = sample(
      x = row.names(x)[ind], 
      size = n_sample, 
      replace = FALSE
    )
    
    message(paste0("        rows choosen to do alignement are ", paste0(list_index[[i_group]], collapse = ",")))
    
    
    # Getting x_M_align: Get the s points from x
    ind_M = which( row.names(x) %in% list_index[[i_group]] )
    sub_matrix =  x[ind_M, ]
    
    if(i_group == 1){
      x_M_align = sub_matrix
    }else{
      x_M_align = rbind(
        x_M_align,
        sub_matrix
      )
    }
  
  }
  
  message(paste0("       After diving in p groups, matrices have length: ", paste0(map_int(list_matrix, nrow), collapse = ",")))
  
  # Up to here we have:
  # x1,..., xp 
  # the list of s poits per each i matrix x in 1:p
  # x_M_align the matrix to align the solutions
  
  # If n/p<l, apply classical MDS
  if( n/p<=l | p == 1){
    message(paste0("Non-recursive call with p: ", p))
    for (j_group in 1:p) {
      distance_matrix = daisy(
        x = list_matrix[[j_group]],
        metric = metric
      )
      
      list_mds[[j_group]] = stats::cmdscale(
        d = distance_matrix, 
        k = s
      )
      
      message(paste0("       j_group is : ", j_group))
      message(paste0("       doing classical MDS with rows : ", paste0( row.names( list_matrix[[j_group]] ), collapse = ',' )))
      message(paste0("       MDS with rows : ", paste0( row.names( list_matrix[[j_group]] ), collapse = ', ')))
      message(paste0("       First component values before alignment are: ", paste0( round( list_mds[[j_group]][, 1], 0 ), collapse = ",") ) )
    }
  }else{
    # Call Fast MDS
    message(paste0("Since n/p>l, calling recursively"))
    
    for(k_group in 1:p){
      message(paste0("       k_group for recursion takes the value ",  k_group))
      message(paste0("       n_rows for recursion ",  nrow(list_matrix[[k_group]])))
      message(paste0("       x matrix for recursion ", paste0(row.names(list_matrix[[k_group]]), collapse = ', ') ))
      fast_mds_2(
        x = list_matrix[[k_group]],
        n = nrow(list_matrix[[k_group]]),
        l = l,
        s = s,
        k = k,
        metric = metric 
      )
    }
  }
  
  total_groups = length(list_mds)
  message(paste0("rows of x_M_align: ", paste0(row.names(x_M_align), collapse = ',')))
  if( total_groups > 0){
    # Build M_align
    # Calculate distance
    distance_matrix_M = daisy(
      x = x_M_align,
      metric = metric
    )
  
    M_align =  stats::cmdscale(
      d = distance_matrix_M, 
      k = s
    ) 
    
    row.names(M_align) = row.names(x_M_align)
  
    # Align the solutions
    message(paste0("There are ", total_groups, " groups to align"))
    
    # Selecting s points points from p submatrices
    for( l_group in 1:total_groups ){
      
      
      # Select the s points from M_align
      rows_names_i = list_index[[l_group]]
      ind_M = which( row.names( M_align ) %in% rows_names_i )
      M_alingn_filter = M_align[ind_M, ]  
      message(paste0("rows from M to align: ", paste0(row.names(M_alingn_filter), collapse = ', ')))
      
    
      # Selecting s points from the MDS of each group
      di = list_mds[[l_group]]
      ind_di = which( row.names(di) %in% list_index[[l_group]])
      di_filter = di[ind_di,]
      message(paste0("rows from D to align: ", paste0(row.names(di_filter), collapse = ', ')))
    
      # Alignment
      procrustes_result =  MCMCpack::procrustes(
        X = di_filter, #The matrix to be transformed
        Xstar = M_alingn_filter, # target matrix
        translation = TRUE, 
        dilation = TRUE
      )
    
      rotation_matrix = procrustes_result$R
      dilation = procrustes_result$s
      translation = procrustes_result$tt
      ones_vector = rep(1, nrow(di)) 
      translation_matrix = ones_vector %*% t(translation)
    
    
      tranformation_di = dilation * di %*% rotation_matrix + translation_matrix
  
      list_Z <<- list.append(
        list_Z,
        tranformation_di
      )
      
      if(is_computed_first_time == TRUE){
        Z_mds <<- tranformation_di
        is_computed_first_time <<- FALSE
      } else{
        Z_mds <<- rbind(
          Z_mds,
          tranformation_di
        )
      }
      
      message(paste0("Value of first row of Z", paste0(Z_mds[1, ], collapse = ",")))
    }
  }
}
   
list_control_index <<- list()
is_computed_first_time <<- TRUE
Z_mds <<- NA
list_row_names <<- list()
list_mds <<- list()
list_Z <<- list()


x = BudgetFood %>% slice(1:100) %>% select(-sex, -town) 
n = nrow(x)
l = 20
s = 2
k = 3
metric = "euclidean"

n*k*s <= l^2

fast_mds_sol <- fast_mds_2(
  x = x,
  n = n,
  l = l,
  s = s,
  k = k,
  metric = "euclidean"
)

nrow(Z_mds)
head(Z_mds)
list_control_index

permutation_order = match(row.names(x), row.names(Z_mds))
Z_mds <<- Z_mds[permutation_order, ]
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



rows choosen to do alignement are 2,7,9,3,6,1
rows choosen to do alignement are 16,17,15,18,20,21


