source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_gower_dist_interpolation <- function(n, l, k) {
  
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

gower_interpolation_formula <- function(A, constant, q_product_all, q_product_last, X_1__S_inv) {
  
  if (nrow(A) == nrow(q_product_all)) {
    q_product <- q_product_all
  } else if (nrow(A) == nrow(q_product_last)) {
    q_product <- q_product_last
  } else {
    stop("Error in dimensions")
  }
  
  new_points <- constant * (q_product - A) %*% X_1__S_inv
  return(new_points)
  
}


extract_rows_cols <- function(matrix, first_rows, last_cols) {
  total_cols <- ncol(matrix)
  return(matrix[1:first_rows, (total_cols - last_cols + 1):total_cols, drop = FALSE])
}

gower_interpolation_dist_mds <- function(x, l, k, dist_fn = stats::dist, ...) {
  
  n <- nrow(x)
  idexes_partition <- get_partitions_for_gower_dist_interpolation(n = n, l = l, k = k)
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
    x_other <- parallel::mclapply(idexes_partition[2:num_partitions],
                                  function(matrix, idx) { matrix[idx, , drop = FALSE] },
                                  matrix = x)
    n_row_x_other <- mclapply(x_other, nrow)
    
    # Obtain the distance matrix with respect the first partition
    m_union <- parallel::mclapply(x_other, function(X, Y) rbind(X, Y), Y = data_1)
    dist_m <- parallel::mclapply(m_union, function(x) as.matrix(dist(x)))

    distance_matrix_filter <- parallel::mcmapply(extract_rows_cols,
                                                 matrix = dist_m,
                                                 first_rows = n_row_x_other,
                                                 MoreArgs = list(last_cols = nrow(data_1)),
                                                 SIMPLIFY = FALSE)
    
    # A matrix
    A <- parallel::mclapply(distance_matrix_filter, function(x){ x^2 })
    ones_vector_all <- matrix(data = 1, nrow = l, ncol = 1)
    ones_vector_last <- matrix(data = 1, nrow = length(idexes_partition[[num_partitions]]), ncol = 1)
    
    
    # Get MDS for all the partitions. The first step is to calculate as many things as possible
    # before performing a lapply
    X_1__S_inv <- X_1 %*% S_inv
    ones_vector_all__q_1_vector <- ones_vector_all %*% t(q_1_vector)
    ones_vector_last__q_1_vector <- ones_vector_last %*% t(q_1_vector)
    constant <- 1/(2*n_1)
    
    MDS <- parallel::mclapply(A, gower_interpolation_formula, 
                              constant = constant,
                              q_product_all = ones_vector_all__q_1_vector,
                              q_product_last = ones_vector_last__q_1_vector,
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
