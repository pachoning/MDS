source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer_setpwise <- function(n, l, tie, k) {
  
  if (l-tie <= 0) {
    stop("l must be greater than tie")
  } else if(l-tie <= tie) {
    stop("l-tie must be greater than tie")
  } else if(l-tie <= k){
    stop("l-tie must be greater than k")
  }
  
  permutation <- sample(x = n, size = n, replace = FALSE)
  
  if (n<=l) {
    list_indexes <- list(permutation)
  } else {
    min_sample_size <- max(k+2, tie)
    p <- 1 + ceiling((n-l)/(l-tie))
    last_partition_sample_size <- n - (l + (p-2) * (l-tie))

    if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
      p <- p - 1
      last_partition_sample_size <- n - (l + (p-2) * (l-tie))
    }

    first_parition <- permutation[1:l]
    last_partition <- permutation[(n-last_partition_sample_size+1):n]
    list_indexes <- split(x = permutation[(l+1):(n-last_partition_sample_size)], f = 1:(p-2))
    names(list_indexes) <- NULL
    list_indexes[[p-1]] <- list_indexes[[1]]
    list_indexes[[p]] <- last_partition
    list_indexes[[1]] <- first_parition
  }

  return(list_indexes)

}


divide_matrix <- function(x, long, idx) {
  if (idx == 1) {
    return(list(first = NA, rest = x))
  } else {
    n_row <- nrow(x)
    x_first <- x[1:long, , drop = FALSE]
    x_rest <- x[(long+1):n_row, , drop = FALSE]
    return(list(first = x_first, rest = x_rest))
  }
}

multidimensional_procrustes <- function(x, target, matrix_to_transform, translation = FALSE, ...) {
  
  if (any(is.na(x)) | any(is.na(target)) | any(is.na(matrix_to_transform))) {
    return(NA)
  } else {
    return(perform_procrustes(x = x, 
                              target = target, 
                              matrix_to_transform = matrix_to_transform,
                              translation = translation,
                              ...))
  }
  
}

divide_conquer_mds <- function(x, l, tie, k, dist_fn = stats::dist, ...) {
  
  n_row_x <- nrow(x)
  
  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    mds_to_return$eigen <- mds_to_return$eigen/nrow(x)
    mds_to_return$GOF <- mds_to_return$GOF
    
  } else {
    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer_setpwise(n = n_row_x, l = l, tie = tie, k = k)
    num_partitions <- length(idx)
    
    # Get partitions
    x_partition <- parallel::mclapply(idx, function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)
    num_rows_partition <- parallel::mclapply(x_partition, nrow)
    
    # Take a sample from the first partition
    idx_sample_1 <- sample(x = num_rows_partition[[1]], size = tie, replace = FALSE)
    x_sample_1 <- x_partition[[1]][idx_sample_1, , drop = FALSE]
    
    # Join each partition with the sample from the first partition except from the first partition
    x_join <- parallel::mclapply(
      x_partition,
      function(m_big, m_small) rbind(m_small, m_big),
      m_small = x_sample_1)
    
    num_rows_x_join <- parallel::mclapply(x_join, nrow)
    x_join[[1]] <- x_join[[1]][(length(idx_sample_1)+1):num_rows_x_join[[1]], , drop = FALSE]
    
    
    # Perform MDS for each partition as well as for the first partition
    mds <- parallel::mclapply(x_join, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE)
    mds_points <- parallel::mclapply(mds, FUN = function(x) x$points) 
    mds_eigen <- parallel::mclapply(mds, FUN = function(x) x$eigen) 
    mds_GOF <- parallel::mclapply(mds, FUN = function(x) x$GOF) 
    
    # For each partition, divide the matrix into two part: 
    # first corresponding to the sample of the first partition
    # the rest of the matrix corresponding to the observations of each matrix
    # Do not divide the first matrix
    mds_division <- parallel::mcmapply(
      divide_matrix,
      x = mds_points,
      idx = 1:num_partitions,
      MoreArgs = list(long = tie),
      SIMPLIFY = FALSE
    )
    
    mds_division_first <- parallel::mclapply(mds_division, function(x) x$first)
    mds_division_rest <- parallel::mclapply(mds_division, function(x) x$rest)
    
    # Get MDS for the sampling poinrs of the first partition
    mds_1_sample <- mds_division_rest[[1]][idx_sample_1, , drop = FALSE]
    
    # Apply Procrustes for each partition
    mds_procrustes <- parallel::mcmapply(FUN = multidimensional_procrustes, 
                                         x = mds_division_first, 
                                         matrix_to_transform = mds_division_rest,
                                         MoreArgs = list(target = mds_1_sample, translation = FALSE),
                                         SIMPLIFY = FALSE)
    
    mds_procrustes[[1]] <- mds_division_rest[[1]]
    
    # Join all the solutions
    mds_solution <- do.call(rbind, mds_procrustes)
    mds_solution <- apply(mds_solution, MARGIN = 2, FUN = function(y) y - mean(y))
    idx_order <- Reduce(c, idx)
    idx_order <- order(idx_order)
    mds_solution <- mds_solution[idx_order, , drop = FALSE]
    mds_solution <- mds_solution %*% eigen(cov(mds_solution))$vectors

    # Get eigenvalues
    eigen <- parallel::mcmapply(function(x, y) x/length(y), x = mds_eigen, y = idx, SIMPLIFY = FALSE)
    eigen <- Reduce(`+`, eigen)
    eigen <- eigen/num_partitions
    
    # Get GOF
    GOF <- parallel::mcmapply(function(x, y) x*length(y), x = mds_GOF, y = idx, SIMPLIFY = FALSE)
    GOF <- Reduce(`+`, GOF)
    GOF <- GOF/n_row_x
    
    mds_to_return <- list(points = mds_solution, eigen = eigen, GOF = GOF)
  }
  
  return(mds_to_return)
}
