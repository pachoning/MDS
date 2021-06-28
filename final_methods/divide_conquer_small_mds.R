source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer_small <- function(n, k, l) {
  
  if (n <= k) {
    stop("n must be greater than k")
  } 

  size_partition <- k + 2
  permutation <- sample(x = n, size = n, replace = FALSE)
  first_sample_size <- l
  
  p <- 1 + ceiling((n-l)/size_partition)
  last_partition_sample_size <- n - (l + (p-2) * size_partition)

  if (last_partition_sample_size < size_partition) {
    p <- p-1
  }
  
  first_partition <- permutation[1:l]
  others_partition <- permutation[(l+1):(l + (p-2)*size_partition)]
  last_partition <- permutation[(l + (p-2)*size_partition + 1):n]
  
  list_indexes <- split(x = others_partition, f = 1:(p-2))
  names(list_indexes) <- NULL
  list_indexes[[p-1]] <- list_indexes[[1]]
  list_indexes[[1]] <- first_partition
  list_indexes[[p]] <- last_partition

  return(list_indexes)
}

divide_data_matrix_small <- function(x, long, idx) {
  if (idx == 1) {
    return(list(first = NA, rest = x))
  } else {
    n_row <- nrow(x)
    x_first <- x[1:long, , drop = FALSE]
    x_rest <- x[(long+1):n_row, , drop = FALSE]
    return(list(first = x_first, rest = x_rest))
  }
}

multidimensional_small_procrustes <- function(x, target, matrix_to_transform, translation = FALSE, ...) {
  
  if (any(is.na(x)) | any(is.na(target)) | any(is.na(matrix_to_transform))) {
    return(NA)
  } else {
    return(perform_procrustes(x = x, 
                              target = target, 
                              matrix_to_transform = matrix_to_transform,
                              translation = translation))
    #return(perform_procrustes(x = x, target = target, matrix_to_transform = matrix_to_transform, translation = translation, ...))
  }
  
}

divide_conquer_small_mds <- function(x, l, tie, k, dist_fn = stats::dist, ...) {

  n_row_x <- nrow(x)
  
  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    #mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen/nrow(x)
    mds_to_return$GOF <- mds_to_return$GOF
    
  } else {

    if(tie > l) {
      stop("l must be greater than tie")
    }

    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer_small(n = n_row_x, k = k, l = l)
    num_partitions <- length(idx)
    
    # Get partitions
    x_partition <- parallel::mclapply(idx, function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)
    num_rows_partition <- parallel::mclapply(x_partition, nrow)
    
    # Take the first partition
    x_1 <- x_partition[[1]]
    
    # Take a sample from the first partition
    idx_1_sample <- sample(x = nrow(x_1), size = tie, replace = FALSE)
    x_1_sample <- x_1[idx_1_sample, , drop = FALSE]
    
    # Join each partition with the first partition except from the first partition
    x_join <- parallel::mclapply(
      x_partition,
      function(m_big, m_small) rbind(m_small, m_big),
      m_small = x_1_sample)
    
    num_rows_x_join <- parallel::mclapply(x_join, nrow)
    x_join[[1]] <- x_join[[1]][(length(idx_1_sample)+1):num_rows_x_join[[1]], , drop = FALSE]
    
    
    # Perform MDS for each partition
    #mds <- parallel::mclapply(x_join, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
    mds <- parallel::mclapply(x_join, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE)
    mds_points <- parallel::mclapply(mds, FUN = function(x) x$points) 
    mds_eigen <- parallel::mclapply(mds, FUN = function(x) x$eigen) 
    mds_GOF <- parallel::mclapply(mds, FUN = function(x) x$GOF) 
    
    # For each partition, divide the matrix into two part: 
    # first corresponding to the first partition
    # the rest of the matrix corresponding to the observations of each matrix
    # Do not divide the first matrix
    mds_division <- parallel::mcmapply(
      divide_data_matrix_small,
      x = mds_points,
      idx = 1:num_partitions,
      MoreArgs = list(long = length(idx_1_sample)),
      SIMPLIFY = FALSE
    )
    
    mds_division_first <- parallel::mclapply(mds_division, function(x) x$first)
    mds_division_rest <- parallel::mclapply(mds_division, function(x) x$rest)
    
    # Get MDS for the first partition
    mds_1_sample <- mds_division_rest[[1]][idx_1_sample, , drop = FALSE]
    
    # Apply Procrustes for each partition
    mds_procrustes <- parallel::mcmapply(FUN = multidimensional_small_procrustes, 
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
    weight <- parallel::mclapply(idx, function(x, sample_size) length(x)/sample_size, sample_size = n_row_x)
    eigen <- parallel::mcmapply(function(x, y) x*y, x = eigen, y = weight, SIMPLIFY = FALSE)
    eigen <- Reduce(`+`, eigen)
    
    # Get GOF
    GOF <- parallel::mcmapply(function(x, y) x*length(y), x = mds_GOF, y = idx, SIMPLIFY = FALSE)
    GOF <- Reduce(`+`, GOF)
    GOF <- GOF/n_row_x
    
    mds_to_return <- list(points = mds_solution, eigen = eigen, GOF = GOF)
  }
  
  return(mds_to_return)
}
