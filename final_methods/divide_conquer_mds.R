source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer <- function(n, l, tie, k) {
  
  size_partition <- l-tie
  p <- ceiling(n/size_partition)
  min_sample_size <- max(k, tie)
  partition_sample_size <- l - tie
  last_sample_size <- n - (p-1)*partition_sample_size
  
  if (l-tie <= 0) {
    stop("l must be greater than tie")
  } else if(l-tie <= tie) {
    stop("l-tie must be greater than tie")
  } else if(l-tie <= k){
    stop("l-tie must be greater than k")
  }
  
  if (last_sample_size < min_sample_size & last_sample_size > 0) {
    p <- p - 1
  }
  
  permutation <- sample(x = p, size = n, replace = TRUE)
  list_indexes <- lapply(1:p, function(x, y) which(x == y), y =  permutation)
  return(list_indexes)

}


divide_matrix <- function(x, long) {
  n_row <- nrow(x)
  x_first <- x[1:long, , drop = FALSE]
  x_rest <- x[(long+1):n_row, , drop = FALSE]
  return(list(first = x_first, rest = x_rest))
}

divide_conquer_mds <- function(x, l, tie, k, dist_fn = stats::dist, ...) {
  
  n_row_x <- nrow(x)
  
  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen/nrow(x)
    mds_to_return$GOF <- mds_to_return$GOF
    
  } else {
    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, tie = tie, k = k)
    num_partitions <- length(idx)
    
    # Get elements of the first partition
    x_1 <- x[idx[[1]], , drop = FALSE]
    
    # Get elements of remaining partitions
    x_rest <- lapply(idx[2:num_partitions], function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)
    
    # Take a sample from the first partition
    idx_sample_1 <- sample(x = nrow(x_1), size = tie, replace = FALSE)
    x_sample_1 <- x_1[idx_sample_1, , drop = FALSE]
    
    # Join each partition with the sample from the first partition
    x_join_1 <- lapply(x_rest, function(m_big, m_small) rbind(m_small, m_big), m_small = x_sample_1)
    
    # Perform MDS for each partition as well as for the first partition
    mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/nrow(x_1)
    mds1_GOF <- mds_1$GOF
    mds_1_sample <- mds_1_points[idx_sample_1, ,drop = FALSE]
    
    mds_join_1 <- lapply(x_join_1, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
    mds_join_1_points <- lapply(mds_join_1, function(x) x$points)
    mds_join_1_eigen <- lapply(mds_join_1, function(x) x$eigen)
    mds_join_1_GOF <- lapply(mds_join_1, function(x) x$GOF)
    
    # For each partition, divide the matrix into two part: 
    # first corresponding to the sample of the first partition
    # the rest of the matrix corresponding to the observations of each matrix
    mds_division <- lapply(mds_join_1_points, divide_matrix, long = tie)
    mds_division_first <- lapply(mds_division, function(x) x$first)
    mds_division_rest <- lapply(mds_division, function(x) x$rest)
    
    # Apply Procrustes for each partition
    mds_procrustes <- mapply(FUN = perform_procrustes, 
                             x = mds_division_first, 
                             matrix_to_transform = mds_division_rest,
                             MoreArgs = list(target = mds_1_sample, translation = FALSE),
                             SIMPLIFY = FALSE)
    
    # Join all the solutions
    mds_solution <- Reduce(rbind, mds_procrustes)
    mds_solution <- Reduce(rbind, list(mds_1_points, mds_solution))
    mds_solution <- apply(mds_solution, MARGIN = 2, FUN = function(y) y - mean(y))
    idx_order <- Reduce(c, idx)
    idx_order <- order(idx_order)
    mds_solution <- mds_solution[idx_order, , drop = FALSE]
    mds_solution <- mds_solution %*% eigen(cov(mds_solution))$vectors

    # Get eigenvalues
    eigen <- mapply(function(x, y) x/length(y), x = mds_join_1_eigen, y = idx[2:num_partitions], SIMPLIFY = FALSE)
    eigen <- Reduce(`+`, eigen)
    eigen <- Reduce(`+`, list(mds_1_eigen, eigen))
    eigen <- eigen/num_partitions
    
    # Get GOF
    GOF_1 <- mds1_GOF*nrow(x_1)
    GOF_rest <- mapply(function(x, y) x*nrow(y), x = mds_join_1_GOF, y = mds_division_rest, SIMPLIFY = FALSE)
    GOF <- Reduce(`+`, GOF_rest)
    GOF <- (GOF + GOF_1)/n_row_x
    
    mds_to_return <- list(points = mds_solution, eigen = eigen, GOF = GOF)
  }
  
  return(mds_to_return)
}
