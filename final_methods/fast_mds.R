source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_fast <- function(n, l, s, k) {
  
  if (n < s) {
    stop("nrow(x) must be greater than s")
  } else if (n*s/l < k) {
    stop("nrow(x)*s/l must be greater than k")
  }
  
  p <- floor(l/s)
  min_sample_size <- max(k+2, s)
  size_partition <- floor(n/p)
  last_sample_size <- n - (p-1) * size_partition
  
  if (last_sample_size < min_sample_size & last_sample_size > 0) {
    p <- p - 1
  }
  
  permutation <- sample(x = n, size = n, replace = FALSE)
  permutation_all <- permutation[1:((p-1)*size_partition)]
  permutation_last <- permutation[((p-1)*size_partition+1):n]
  list_indexes <- split(x = permutation_all, f = 1:(p-1))
  names(list_indexes) <- NULL
  list_indexes[[p]] <- permutation_last
  
  return(list_indexes)
}

#main_fast_mds <- function(idx, matrix, l, s, k, dist_fn, n_cores, ...) {
main_fast_mds <- function(idx, matrix, l, s, k, dist_fn, n_cores) {
  
  # Partition the matrix
  x_partition <- matrix[idx, , drop = FALSE]
  
  # Apply the method
  #mds <- fast_mds(x = x_partition, l = l, s = s, k = k, dist_fn = dist_fn, n_cores = n_cores, ...)
  mds <- fast_mds(x = x_partition, l = l, s = s, k = k, dist_fn = dist_fn, n_cores = n_cores)
  
  return(mds)
  
}

fast_mds <- function(x, l, s, k, dist_fn = stats::dist, n_cores = 1, ...) {

  n <- nrow(x)

  if (n <= l) {

    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    #mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds$eigen <- mds$eigen / nrow(x)
    return(mds)

  } else {

    # Split x
    index_partition <- get_partitions_for_fast(n = n, l = l, s = s, k = k)
    num_partition <- length(index_partition)
                                      
    # Apply MDS to all the partitions
#    mds_partition <- parallel::mclapply(index_partition, 
#                                        main_fast_mds,
#                                        matrix = x,
#                                        l = l, 
#                                        s = s, 
#                                        k = k, 
#                                        dist_fn = dist_fn, n_cores = n_cores,
#                                        mc.cores = n_cores,
#                                        ...)

    # Get MDS for each partition recursevely
    mds_partition <- parallel::mclapply(index_partition, 
                                        main_fast_mds,
                                        matrix = x,
                                        l = l, 
                                        s = s, 
                                        k = k, 
                                        dist_fn = dist_fn,
                                        n_cores = n_cores,
                                        mc.cores = n_cores)

    mds_partition_points <- parallel::mclapply(mds_partition, function(x) x$points, mc.cores = n_cores)
    mds_partition_eigen <- parallel::mclapply(mds_partition, function(x) x$eigen, mc.cores = n_cores)
    mds_GOF <- parallel::mclapply(mds_partition, function(x) x$GOF, mc.cores = n_cores)
    
    # take a sample for each partition
    length_partition <- parallel::mclapply(index_partition, length, mc.cores = n_cores)
    sample_partition <- parallel::mclapply(length_partition, sample, size = s, replace = FALSE, mc.cores = n_cores)
    indexes_filtered <- parallel::mcmapply(function(idx, sample) idx[sample],
                                           idx = index_partition,
                                           sample = sample_partition,
                                           SIMPLIFY = FALSE,
                                           mc.cores = n_cores)
    
    length_sample <- parallel::mclapply(sample_partition, length, mc.cores = n_cores)
    
    indexes_scaled <- parallel::mcmapply(function(i, long) ((i-1)*long + 1):(i*long),
                                         i = 1:num_partition,
                                         long = length_sample,
                                         SIMPLIFY = FALSE,
                                         mc.cores = n_cores)

    # Join all the points
    x_partition_sample <- parallel::mclapply(indexes_filtered,
                                             function(index_partitions, matrix) { matrix[index_partitions, , drop = FALSE] },
                                             matrix = x,
                                             mc.cores = n_cores)

    x_M <- do.call(rbind, x_partition_sample)

    # Apply MDS to the subsampling points
    #mds_M <- classical_mds(x = x_M, k = k, dist_fn = dist_fn, ...)
    mds_M <- classical_mds(x = x_M, k = k, dist_fn = dist_fn)
    mds_M_points <- mds_M$points

    # Extract the MDS configuration for the sampling points from mds_M_points 
    mds_M_sampling_points <- parallel::mclapply(indexes_scaled,
                                                function(indexes_scaled, matrix) { matrix[indexes_scaled, , drop = FALSE] },
                                                matrix = mds_M_points,
                                                mc.cores = n_cores)

    # Extract the MDS configuration for the sampling points from mds_partition_points
    mds_partition_sampling_points <- parallel::mcmapply(function(matrix, index_partitions, idx) { matrix[idx, , drop = FALSE] },
                                                        matrix = mds_partition_points,
                                                        idx = sample_partition,
                                                        SIMPLIFY = FALSE,
                                                        mc.cores = n_cores)

    # Apply Procrustes
    procrustes <- parallel::mcmapply(perform_procrustes, 
                                     x = mds_partition_sampling_points, 
                                     target = mds_M_sampling_points, 
                                     matrix_to_transform = mds_partition_points, 
                                     translation = FALSE, 
                                     SIMPLIFY = FALSE,
                                     mc.cores = n_cores)
    
    # Build the list to be returned
    idx_order <- Reduce(c, index_partition)
    idx_order <- order(idx_order)
    mds <-do.call(rbind, procrustes)
    mds <- mds[idx_order, ,drop = FALSE]
    mds <- apply(mds, MARGIN = 2, FUN = function(y) y - mean(y))
    mds <- mds %*% eigen(cov(mds))$vectors
    eigen <- Reduce(`+`, mds_partition_eigen)/num_partition
    
    # Build GOF metric
    GOF <- parallel::mcmapply(function(x, y) x*length(y), 
                              x = mds_GOF, 
                              y = index_partition, 
                              SIMPLIFY = FALSE, 
                              mc.cores = n_cores)
    GOF <- Reduce(`+`, GOF)/n
    
    return(list(points = mds, eigen = eigen, GOF = GOF))
  }
}

