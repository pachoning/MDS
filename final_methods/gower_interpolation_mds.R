source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_gower_interpolation <- function(n, l, k) {
  
  if (l<=k) {
    stop("l must be greater than k")
  }
  p <- ceiling(n/l)
  p <- pmax(1, p)
  
  permutation <- sample(x = n, size = n, replace = FALSE)
  
  if (p>1) {
    p <- ceiling(n/l)
    fix_elements <- permutation[1:(l*(p-1))]
    residual_elements <- permutation[(l*(p-1) + 1):n]
    list_index <- split(x = fix_elements, f = 1:(p-1))
    names(list_index) <- NULL
    list_index[[p]] <- residual_elements

  } else {
    list_index <- list(permutation)
  }
  
  return(list_index)
}


gower_interpolation_mds <- function(x, l, k, dist_fn = stats::dist, ...) {
  
  n <- nrow(x)
  idexes_partition <- get_partitions_for_gower_interpolation(n = n, l = l, k = k)
  num_partitions <- length(idexes_partition)
  
  
  if (num_partitions <= 1) {
    
    # It is possible to run MDS directly
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    points <- mds$points
    eigen <- mds$eigen/n
    GOF <- mds$GOF
    list_to_return <- list(points = points, eigen = eigen, GOF = GOF)
    
  } else {
    
    # Get the first group 
    n_1 <- length(idexes_partition[[1]])

    # Obtain MDS for the first group
    x_1 <- x[idexes_partition[[1]], ,drop = FALSE]
    mds_eig <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    distance_matrix <- mds_eig$distance
    
    M <- mds_eig$points
    eigen <- mds_eig$eigen / nrow(M)
    GOF <- mds_eig$GOF
    
    # Calculations needed to do Gower interpolation
    delta_matrix <- distance_matrix^2
    I_l <- diag(n_1)
    ones_vector <- matrix(data = 1, nrow = n_1, ncol = 1)
    P <- I_l - 1 / n_1 * ones_vector %*% t(ones_vector)
    Q_1 <- -1 / 2 * P %*% delta_matrix %*% t(P) 
    q_1_vector <- diag(Q_1)
    S <- 1 / (nrow(M)-1) * t(M) %*% M
    S_inv <- solve(S)
    
    # Get x for each partition
    x_other <- lapply(idexes_partition[2:num_partitions],
                      function(matrix, idx) { matrix[idx, , drop = FALSE] },
                      matrix = x)
    
    # Obtain the distance matrix with respect the first partition
    distance_matrix_filter <- lapply(x_other, function(X, Y){ pdist::pdist(X, Y) }, Y = x_1)
    distance_matrix_filter <- lapply(distance_matrix_filter, as.matrix)
    
    # A matrix
    A <- lapply(distance_matrix_filter, function(x){ x^2 })
    ones_vector <- lapply(idexes_partition[2:num_partitions], 
                          function(times, x){ rep(x, length(times)) }, 
                          x = 1)
    
    # Get MDS for all the partitions
    MDS <- mapply(function(A, ones_vector) { 1 / (2 * n_1) * (ones_vector %*% t(q_1_vector) - A) %*% M %*% S_inv  }, 
                  A = A, ones_vector = ones_vector, SIMPLIFY = FALSE)
    
    # Get cummulative MDS
    cum_mds <- Reduce(rbind, MDS)
    cum_mds <- Reduce(rbind, list(M, cum_mds))
    
    # Reorder the rows
    idexes_order <- Reduce(c, idexes_partition)
    cum_mds <- cum_mds[order(idexes_order), , drop = FALSE]
    cum_mds <- apply(cum_mds, MARGIN = 2, FUN = function(y) y - mean(y)) 
    cum_mds <- cum_mds %*% eigen(cov(cum_mds))$vectors
    
    # List to return
    list_to_return <- list(points = cum_mds, eigen = eigen, GOF = GOF)
  }
  
  return(list_to_return)
}
