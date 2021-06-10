source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_gower_random_partitions <- function(n, l) {

  # We first choose l indices at random
  idx_1 <- sample(x = n, size = l, replace = FALSE)
  
  
  # Get rest of indexes
  idx_rest <- setdiff(1:n, idx_1)
  
  return(list(first = idx_1, others = idx_rest))
  
}


get_gower_coordinates <- function(A_sq, q_1, S_1_inv, X_1, l){
  m <- dim(A_sq)[1]
  one_vector <- matrix(data = 1, nrow = m, ncol = 1)
  X_2 <- 1/(2*l)*(one_vector %*% t(q_1) - A_sq) %*% X_1 %*% S_1_inv
  return(X_2)
}

gower_by_row <- function(x, l, k, dist_fn = stats::dist, ...) {
  
  n <- nrow(x)
  
  if (n<=l) {
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
  } else {
    # Get indexes partitions
    idx_partitions <- get_gower_random_partitions(n = n, l = l)
    idx_first <- idx_partitions$first
    idx_others <- idx_partitions$others
    idx_others <- as.list(idx_others)
    
    # Get the matrix for the first partition
    x_1 <- x[idx_first, , drop = FALSE]
    
    # Get the set of matrices
    x_others <- lapply(idx_others, function(indexes, matrix) matrix[indexes, , drop = FALSE], matrix = x)
    
    
    # Perform mds over the first partition
    mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    points_1 <- mds_1$points
    eigen_1 <- mds_1$eigen/nrow(points_1)
    GOF_1 <- mds_1$GOF
    D_dist <- mds_1$distance
    
    # Get the matrix A for all the points
    A <- lapply(x_others, function(matrix_1, matrix_2) pdist::pdist(matrix_1, matrix_2), matrix_2 = x_1)
    A <- lapply(A, as.matrix)
    
    # Get Delta_1 matrix
    Delta_1 <- D_dist^2
    
    # Get P matrix
    P <- diag(l) - 1/l * matrix(data=1, nrow = l, ncol = 1) %*% matrix(data = 1, nrow = 1, ncol = l)
    
    # Get Q_1 matrix
    Q_1 <- -1/2*P %*% Delta_1 %*% t(P)
    q_1 <- matrix(data = diag(Q_1), nrow = nrow(Q_1), ncol = 1)
    
    # Get sqaured of A
    A_sq <- lapply(A, function(x) x^2)
    
    # Get S_1
    S_1 <- 1/(nrow(points_1) - 1) * t(points_1) %*% points_1
    S_1_inv <- solve(S_1)
    
    # Get configurations for the remaining parts
    points_others <- mapply(get_gower_coordinates, 
                            A_sq, 
                            MoreArgs = list(q_1 = q_1, S_1_inv = S_1_inv, X_1 = points_1, l = l),
                            SIMPLIFY = FALSE)
    
    # Append all configurations
    mds_append <- Reduce(rbind, points_others)
    mds_append <- rbind(points_1, mds_append)
    
    # Rearrange indexes
    idx_all <- c(idx_partitions$first, idx_partitions$others)
    idx_order <- order(idx_all)
    mds_config <- mds_append[idx_order, , drop = FALSE]
    mds_config <- apply(mds_config, MARGIN = 2, FUN = function(y) y - mean(y)) 
    mds_config <- mds_config %*% eigen(cov(mds_config))$vectors
    mds <- list(points = mds_config, eigen = eigen_1, GOF = GOF_1)
  }
  
  
  return(mds)
}
