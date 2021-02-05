source("tools/classical_mds.R")
source("tools/procrustes.R")


divide_conquer_single_partiton <- function(x_other, x_first_sample, mds_first, k, dist_fn, ...) {
  size_first_sample <- nrow(x_first_sample)
  x_join <- rbind(x_first_sample, x_other)
  n_row_x_join <- nrow(x_join)
  mds_x_join <- classical_mds(x = x_join, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
  mds_x_join_points <- mds_x_join$points
  mds_x_join_eigen <- mds_x_join$eigen/n_row_x_join
  mds_x_join_first <- mds_x_join_points[1:3, , drop = FALSE]
  mds_x_join_other <- mds_x_join_points[(size_first_sample + 1):n_row_x_join, , drop = FALSE]
  mds_other_aligned <- perform_procrustes(x = mds_x_join_first, target = mds_first, matrix_to_transform = mds_x_join_other, 
                                          translation = FALSE, dilation = FALSE)
  
  return(list(points = mds_other_aligned, eigen = mds_x_join_eigen))
}


divide_conquer_mds_efficient <- function(x, l, tie, k, dist_fn = stats::dist, ...) {

  n_row_x <- nrow(x)

  if (n_row_x <= l) {

    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen / nrow(x)
    
  } else {
    # Generate indexes list. Each element correspond to the index of the partition
    init_time <- proc.time()
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, tie = tie, k = k)
    elapsed_time <- proc.time() - init_time
    message(paste0("  Time to create the partition: ", elapsed_time[3]))
    num_partitions <- length(idx)
    
    # Get the elements of the first partition and drop from the list
    idx_first <- idx[[1]]
    idx[[1]] <- NULL 
    
    # Partition x
    x_partition <- lapply(idx, function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)
    
    # Take first partition
    x_first <- x[idx_first, , drop = FALSE]
    
    # Take a sample from the first partition
    sample_first_partition <- sample(x = idx_first, size = tie, replace = FALSE)
    x_first_sample <- x_first[sample_first_partition, , drop = FALSE]
    
    # Perform MDS for the first partition
    mds_first <- classical_mds(x = x_first, k = k, dist_fn = dist_fn)
    mds_first_points <- mds_first$points
    mds_first_eigen <- mds_first$eigen/nrow(x_first)
    
    
    # Obtain MDS for the remaining partitions
    mds_others <- lapply(x_partition, divide_conquer_single_partiton, 
                         x_first_sample = x_first_sample, mds_first = mds_first_points, k = k, dist_fn = dist_fn, ...)
    
    mds_others_points <- lapply(mds_others, function(x) x$points)
    mds_others_eigen <- lapply(mds_others, function(x) x$eigen)
    
    # Build MDS for the entire matrix
    mds_solution <- Reduce(rbind, mds_others_points)
    mds_solution <- Reduce(rbind, list(mds_first_points, mds_solution))
    
    # Build eigenvalues for the entire matrix
    eigen <- Reduce(`+`, list(mds_first_eigen, Reduce(`+`, mds_others_eigen)))
    eigen <- eigen/num_partitions
    
    # Build the list to be returned
    mds_to_return <- list(points = mds_solution, eigen = eigen)
  }
  return(mds_to_return)
}