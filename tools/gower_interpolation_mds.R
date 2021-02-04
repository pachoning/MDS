source("tools/classical_mds.R")
source("tools/procrustes.R")

gower_interpolation_mds <- function(x, l, k, dist_fn = stats::dist, ...) {
  
  nrow_x <- nrow(x)
  p <- ceiling(nrow_x / l)
  
  if (p < 1) {
    p <- 1
  } 
  
  if (p > 1) {
    
    initial_row_names <- row.names(x)
    row.names(x) <- 1:nrow(x)
    
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution <- sample(x = p, size = nrow_x, replace = TRUE)
    
    # Get the first group 
    ind_1 <- which(sample_distribution == 1)
    n_1 <- length(ind_1)
    
    # Do MDS with the first group
    submatrix_data <- x[ind_1, ,drop = FALSE]
    mds_eig <- classical_mds(x = submatrix_data, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    distance_matrix <- mds_eig$distance
    
    M <- mds_eig$points
    eigen <- mds_eig$eigen / nrow(M)
    cum_mds <- M
    
    # Calculations needed to do Gower interpolation
    delta_matrix <- distance_matrix^2
    In <- diag(n_1)
    ones_vector <- rep(1, n_1)
    J <- In - 1 / n_1 * ones_vector %*% t(ones_vector)
    G <- -1 / 2 * J %*% delta_matrix %*% t(J) 
    g_vector <- diag(G)
    S <- 1 / (nrow(M)-1) * t(M) %*% M
    S_inv <- solve(S)
    
    # For the rest of the groups, do the interpolation
    for (i_group in 2:p) {
      # Filtering the data
      ind_i_group <- which(sample_distribution == i_group)
      submatrix_data <- x[ind_i_group, ,drop = FALSE]
      
      # A matrix
      distance_matrix_filter <- pdist::pdist(
        X = submatrix_data,
        Y = x[ind_1, ,drop = FALSE]
      )
      
      distance_matrix_filter <- as.matrix(distance_matrix_filter)
      A <- distance_matrix_filter^2
      ones_vector <- rep(1, length(ind_i_group))
      MDS_i_group <- 1 / (2 * n_1) * (ones_vector %*% t(g_vector) - A) %*% M %*% S_inv
      row.names(MDS_i_group) <- row.names(submatrix_data)
      cum_mds <- rbind(cum_mds, MDS_i_group)
    }
    
    cum_mds <- cum_mds[order(as.numeric(row.names(cum_mds))), , drop = FALSE]
    row.names(cum_mds) <- initial_row_names
    row.names(x) <- initial_row_names
    
  } else {
    # It is possible to run MDS directly
    mds_eig <- classical_mds(x = x, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    distance_matrix <- mds_eig$distance
    cum_mds <- mds_eig$points
    eigen <- mds_eig$eigen / nrow_x
  }
  
  return(list(points = cum_mds, eigen = eigen))
}
