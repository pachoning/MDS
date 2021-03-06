source("tools/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_fast <- function(n, l, s, k) {

  if (n/l<1.1) {
    p <- 2
  } else {
    p <- ceiling(l/s)
  }

  min_sample_size <- max(k+2, s)
  size_partition <- floor(n/p)
  last_sample_size <- n - (p-1)*size_partition
  list_indexes <- list()
  
  if (size_partition < s) {
    stop("nrow(x) must be greater than s")
  } else if (size_partition < k) {
    stop("nrow(x)*s/l must be greater than k")
  }
  
  if (last_sample_size < min_sample_size & last_sample_size > 0) {
    p <- p - 1
  }
  
  for(i in 1:p) {
    if (i == 1) {
      ini <- 1
      end <- size_partition
    } else if (i < p) {
      ini <- end + 1
      end <- (ini - 1) + size_partition
    } else {
      ini <- end + 1
      end <- n
    }
    
    list_indexes[[i]] <- ini:end
  }
  
  return(list_indexes)
}



fast_mds <- function(x, l, s, k, dist_fn = stats::dist, ...) {
  
  n <- nrow(x)
  
  if (n <= l) {
    
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds$eigen <- mds$eigen / nrow(x)
    return(mds)
    
  } else {
    
    # Split x
    index_partition <- get_partitions_for_fast(n = n, l = l, s = s, k = k)
    x_partition <- lapply(index_partition, function(idx, matrix) { matrix[idx, , drop = FALSE] }, matrix = x)
    num_partition <- length(index_partition)
    
    # Apply MDS to all the partitions
    mds_partition <- lapply(x_partition, fast_mds, l = l, s = s, k = k, dist_fn = dist_fn, ...)
    mds_partition_points <- lapply(mds_partition, function(x) x$points)
    mds_partition_eigen <- lapply(mds_partition, function(x) x$eigen)
    
    # take a subsample for each partition
    length_partition <- lapply(index_partition, length)
    sample_partition <- lapply(length_partition, sample, size = s, replace = FALSE)
    x_partition_sample <- mapply(function(matrix, idx) { matrix[idx, , drop = FALSE] }, 
                                 matrix = x_partition, idx = sample_partition, SIMPLIFY = FALSE)
    
    
    # Create two lists: one with initial position inside the matrix and another with end position
    ini_index <- list()
    end_index <- list()
    
    for (i in 1:num_partition) {
      length_sample_partition_i <- length(sample_partition[[i]])
      ini_index[[i]] <- (i-1)*length_sample_partition_i + 1
      end_index[[i]] <- i*length_sample_partition_i
    }
    
    # Join each sampled data
    x_M <- Reduce(rbind, x_partition_sample)
    
    # Apply MDS to the subsampling points
    mds_M <- classical_mds(x = x_M, k = k, dist_fn = stats::dist, ...)
    mds_M_points <- mds_M$points
    
    # Extract the MDS configuration for the sampling points from mds_M_points 
    mds_M_sampling_points <- mapply(function(matrix, ini, end) {  matrix[ini:end, , drop = FALSE] }, 
                                    ini = ini_index, end = end_index, 
                                    MoreArgs = list(matrix = mds_M_points), SIMPLIFY = FALSE)
    
    # Extract the MDS configuration for the sampling points from mds_partition_points
    mds_partition_sampling_points <- mapply(function(matrix, idx) { matrix[idx, , drop = FALSE] }, 
                                            matrix = mds_partition_points, idx = sample_partition, SIMPLIFY = FALSE)
    
    # Apply Procrustes
    procrustes <- mapply(perform_procrustes, x = mds_partition_sampling_points, target = mds_M_sampling_points, 
                         matrix_to_transform = mds_partition_points, translation = FALSE, dilation = FALSE, SIMPLIFY = FALSE)
    
    # Build the list to be returned
    mds <- Reduce(rbind, procrustes)
    eigen <- Reduce(`+`, mds_partition_eigen)/num_partition
    
    return(list(points = mds, eigen = eigen))
  }
}
