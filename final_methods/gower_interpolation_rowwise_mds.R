source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_gower_interpolation_rowwise <- function(n, l) {
  
  if(n <= l) {
    p <- 1
  } else {
    p <- 1 + (n - l)
  }
  
  permutation <- sample(x = n, size = n, replace = FALSE)
  
  if (p>1) {
    first_part <- permutation[1:l]
    remaining_part <- permutation[(l+1):n]
    list_index <- split(x = remaining_part, f = 1:(p-1))
    names(list_index) <- NULL
    list_index[[p]] <- list_index[[1]]
    list_index[[1]] <- first_part
  } else {
    list_index <- list(permutation)
  }

  return(list_index)
}

gower_interpolation_formula_rowwise <- function(A, constant, q_product, X_1__S_inv) {
  
  new_points <- constant * (q_product - A) %*% X_1__S_inv
  return(new_points)
  
}


gower_interpolation_rowwise_mds <- function(x, l, k, dist_fn = stats::dist, ...) {
  
  n <- nrow(x)
  idexes_partition <- get_partitions_for_gower_interpolation_rowwise(n = n, l = l)
  num_partitions <- length(idexes_partition)
  
  
  if (num_partitions <= 1) {
    
    # It is possible to run MDS directly
    #mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    points <- mds$points
    eigen <- mds$eigen/n
    GOF <- mds$GOF
    list_to_return <- list(points = points, eigen = eigen, GOF = GOF)
    
  } else {
    
    # Get the first group 
    n_1 <- length(idexes_partition[[1]])
    
    # Obtain MDS for the first group
    data_1 <- x[idexes_partition[[1]], ,drop = FALSE]
    #mds_eig <- classical_mds(x = data_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    mds_eig <- classical_mds(x = data_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE)
    distance_matrix <- mds_eig$distance
    
    X_1 <- mds_eig$points
    eigen <- mds_eig$eigen / nrow(X_1)
    GOF <- mds_eig$GOF
    
    # Calculations needed to do Gower interpolation
    delta_matrix <- distance_matrix^2
    I_l <- diag(n_1)
    ones_vector <- matrix(data = 1, nrow = n_1, ncol = 1)
    P <- I_l - 1 / n_1 * ones_vector %*% t(ones_vector)
    Q_1 <- -1 / 2 * P %*% delta_matrix %*% t(P) 
    q_1_vector <- diag(Q_1)
    S <- 1 / (nrow(X_1)-1) * t(X_1) %*% X_1
    S_inv <- solve(S)
    
    # Get x for each partition
    x_other <- lapply(idexes_partition[2:num_partitions],
                      function(matrix, idx) { matrix[idx, , drop = FALSE] },
                      matrix = x)
    
    # Obtain the distance matrix with respect the first partition
    distance_matrix_filter <- parallel::mclapply(x_other, function(X, Y){ pdist::pdist(X, Y) }, Y = data_1)
    distance_matrix_filter <- parallel::mclapply(distance_matrix_filter, as.matrix)
    
    # A matrix
    A <- parallel::mclapply(distance_matrix_filter, function(x){ x^2 })
    ones_vector_all <- matrix(data = 1, nrow = 1, ncol = 1)

    # Get MDS for all the partitions. The first step is to calculate as many things as possible
    # before performing a lapply
    X_1__S_inv <- X_1 %*% S_inv
    ones_vector_all__q_1_vector <- ones_vector_all %*% t(q_1_vector)
    constant <- 1/(2*n_1)
    
    MDS <- parallel::mclapply(A, 
                  gower_interpolation_formula_rowwise, 
                  constant = constant,
                  q_product = ones_vector_all__q_1_vector,
                  X_1__S_inv = X_1__S_inv)
    
    # Get cummulative MDS
    #cum_mds <- Reduce(rbind, MDS)
    #cum_mds <- Reduce(rbind, list(X_1, cum_mds))
    cum_mds <- do.call(rbind, MDS)
    cum_mds <- rbind(X_1, cum_mds)
    
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
