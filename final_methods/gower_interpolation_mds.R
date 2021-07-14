source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_gower_interpolation <- function(n, n_obs, l, k) {

  if (l<=k) {
    stop("l must be greater than k")
  }

  p <- 1 + ceiling((n - l)/n_obs)
  n_last <- n - (l + (p-2)*n_obs)

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (p>1) {
    first_part <- permutation[1:l]
    middle_part <- permutation[(l+1):(n-n_last)]
    last_part <- permutation[(n-n_last+1):n]
    
    list_index <- split(middle_part, 1:(p-2))
    names(list_index) <- NULL
    list_index[[p-1]] <- list_index[[1]]
    list_index[[1]] <- first_part
    list_index[[p]] <- last_part

  } else {
    list_index <- list(permutation)
  }
  
  return(list_index)
}

get_P_matrix <- function(n_row) {
  
  identity_matrix <- diag(x = 1, nrow = n_row, ncol = n_row)
  one_vector <- matrix(data = 1, nrow = n_row, ncol = 1)
  P <- identity_matrix - 1/n_row * one_vector %*% t(one_vector)
  return(P)

}

gower_mds_main <- function(idx, x, data_1, x_1, n_row_1, q_vector, x_1__s_1__inv) {

  # Filter the matrix
  x_other <- x[idx, , drop = FALSE]
  n_row_other <- nrow(x_other)
  
  # Get A matrix
  A <- as.matrix(pdist(x_other, data_1))
  
  # Get delta matrix
  A_sq <- A^2
  
  # One vecotr
  one_vector_other <- matrix(data = 1, nrow = n_row_other, ncol = 1)
  
  # Get coordinates
  x_2 <- 1/(2*n_row_1) * (one_vector_other %*% t(q_vector) - A_sq) %*% x_1__s_1__inv
  
  return(x_2)
}

gower_interpolation_mds <- function(x, l, k, dist_fn = stats::dist, n_row_partition = l, n_cores = 1, ...) {

  n <- nrow(x)
  indexes_partition <- get_partitions_for_gower_interpolation(n = n, n_obs = n_row_partition, l = l, k = k)
  num_partitions <- length(indexes_partition)
  
  
  if (num_partitions <= 1) {
    
    # It is possible to run MDS directly
    #mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    points <- mds$points
    eigen_v <- mds$eigen/n
    GOF <- mds$GOF
    list_to_return <- list(points = points, eigen = eigen_v, GOF = GOF)
    
  } else {
    
    # Get the first group 
    n_row_1 <- length(indexes_partition[[1]])

    # Obtain MDS for the first group
    data_1 <- x[indexes_partition[[1]], ,drop = FALSE]
    #mds_eig <- classical_mds(x = data_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    mds_eig <- classical_mds(x = data_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE)
    distance_matrix <- as.matrix(mds_eig$distance)
    X_1 <- mds_eig$points
    eigen_v <- mds_eig$eigen/nrow(X_1)
    GOF <- mds_eig$GOF
    
    # Get P matrix
    P <- get_P_matrix(n_row = n_row_1)
    Q <- -1/2 * P %*% distance_matrix^2 %*% t(P)
    q_vector <- diag(Q)
    S <- 1 / (n_row_1-1) * t(X_1) %*% X_1
    x_1__s_1__inv <- X_1 %*% solve(S)

    # Calculations needed to do Gower interpolation
    mds_others <- parallel::mclapply(indexes_partition[2:num_partitions],
                                     gower_mds_main, 
                                     x = x, 
                                     data_1 = data_1,
                                     x_1 = X_1, 
                                     n_row_1 = n_row_1, 
                                     q_vector = q_vector,
                                     x_1__s_1__inv = x_1__s_1__inv,
                                     mc.cores = n_cores)
    
    mds_points <- do.call(rbind, mds_others)
    mds_points <- rbind(X_1, mds_points)
    idx_all <- do.call(c, indexes_partition)
    idx_all <- order(idx_all)
    mds_points <- mds_points[idx_all, , drop = FALSE]
    mds_points <- apply(mds_points, MARGIN = 2, FUN = function(y) y - mean(y))
    mds_points <- mds_points %*% eigen(cov(mds_points))$vectors
    list_to_return <- list(points = mds_points, eigen = eigen_v, GOF = GOF)
  }
  
  return(list_to_return)
}
