gower.interpolation.mds <- function(
  x,
  l,
  s
){
  
  nrow_x = nrow(x)
  p = ceiling(nrow_x/l)
  if(p<1) p = 1
  
  if( p>1 ){
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution = sample(x = p, size = nrow_x, replace = TRUE)
    sample_distribution = sort(sample_distribution)
    
    # Get the first group 
    ind_1 = which(sample_distribution == 1)
    n_1 = length(ind_1)
    
    # Do MDS with the first group
    submatrix_data = x[ind_1, ]
    
    distance_matrix = cluster::daisy(
      x = submatrix_data,
      metric = "euclidean"
    )
    
    distance_matrix = as.matrix(distance_matrix)
    
    # MDS for the first group
    cmd_eig = stats::cmdscale(
      d = distance_matrix, 
      k = s,
      eig = TRUE
    )
    
    M = cmd_eig$points
    eigen = cmd_eig$eig/nrow(M)
    cum_mds = M
    
    # Calculations needed to do Gower interpolation
    delta_matrix = distance_matrix^2 
    In = diag(n_1)
    ones_vector = rep(1, n_1)
    J = In - 1/n_1*ones_vector %*% t(ones_vector)
    G = -1/2 * J %*% delta_matrix %*% t(J) 
    g_vector = diag(G)
    # S = cov(M)
    S = 1/(nrow(M)-1)*t(M) %*% M
    S_inv = solve(S)
    
    # For the rest of the groups, do the interpolation
    for(i_group in 2:p){
      # Filtering the data
      ind_i_group = which(sample_distribution == i_group)
      submatrix_data = x[ind_i_group, ]
      
      
      # A matrix
      distance_matrix_filter = pdist::pdist(
        X = submatrix_data,
        Y = x[ind_1, ]
      )
      
      distance_matrix_filter = as.matrix(distance_matrix_filter)
      A = distance_matrix_filter^2
      ones_vector = rep(1, length(ind_i_group))
      MDS_i_group = 1/(2*n_1)*(ones_vector %*%t(g_vector) - A) %*% M %*% S_inv
      cum_mds = rbind(
        cum_mds,
        MDS_i_group
      )
    }
  }else{
    # It is possible to run MDS directly
    distance_matrix = cluster::daisy(
      x = x,
      metric = "euclidean"
    )
    
    distance_matrix = as.matrix(distance_matrix)
    
    # MDS for the first groups
    cmd_eig = stats::cmdscale(
      d = distance_matrix, 
      k = s,
      eig = TRUE
    )
    
    cum_mds = cmd_eig$points
    eigen = cmd_eig$eig/nrow_x
  }
  
  return(
    list(
      points = cum_mds,
      eig = eigen
    )
  )
}
