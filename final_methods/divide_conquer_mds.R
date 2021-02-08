source("tools/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer <- function(n, l, tie, k) {
  
  size_partition <- l-tie
  p <- ceiling(n/size_partition)
  min_sample_size <- max(k, tie)
  partition_sample_size <- l - tie
  last_sample_size <- n - (p-1)*partition_sample_size
  list_indexes <- list()
  
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
  
  for(i in 1:p){
    if (i == 1) {
      ini <- 1
      end <- size_partition
    } else if (i < p) {
      ini <- end + 1
      end <- (ini-1) + size_partition
    } else {
      ini <- end + 1
      end <- n
    }
    
    list_indexes[[i]] <- ini:end
  }
  
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
    mds_to_return$eigen <- mds_to_return$eigen / nrow(x)
    
  } else {
    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, tie = tie, k = k)
    num_partitions <- length(idx)
    
    # Get the elements of the first partition and drop from the list
    idx_1 <- idx[[1]]
    idx[[1]] <- NULL 
    
    # Partition x
    x_partition <- lapply(idx, function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)
    
    # Take a sample from the first partition
    sample_first_partition <- sample(x = idx_1, size = tie, replace = FALSE)
    x_1 <- x[idx_1, , drop = FALSE]
    x_sample_1 <- x_1[sample_first_partition, , drop = FALSE]
    
    # Join each partition with the sample from the first partition
    x_join_1 <- lapply(x_partition, function(m_big, m_small) rbind(m_small, m_big), m_small = x_sample_1)
    
    # Perform MDS for each partition as well as for the first partition
    mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/nrow(x_1)
    mds_1_sample <- mds_1_points[sample_first_partition, ,drop = FALSE]
    
    mds_join_1 <- lapply(x_join_1, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE)
    mds_join_1_points <- lapply(mds_join_1, function(x) x$points)
    mds_join_1_eigen <- lapply(mds_join_1, function(x) x$eigen)
    
    
    # For each partition, divide the matrix into two part: 
    ## first corresponding to the sample of the first partition
    # the rest of the matrix corresponding to the observations of each matrix
    mds_division <- lapply(mds_join_1_points, divide_matrix, long = tie)
    mds_division_first <- lapply(mds_division, function(x) x$first)
    mds_division_rest <- lapply(mds_division, function(x) x$rest)
    
    # Apply Procrustes for each partition
    mds_procrustes <- mapply(FUN = perform_procrustes, 
                             x = mds_division_first, 
                             matrix_to_transform = mds_division_rest,
                             MoreArgs = list(target = mds_1_sample, translation = FALSE, dilation = FALSE))
    
    # Join all the solutions
    mds_solution <- Reduce(rbind, mds_procrustes)
    mds_solution <- Reduce(rbind, list(mds_1_points, mds_solution))
    
    # Get eigenvalues
    eigen <- mapply(function(x, y) x/length(y), x = mds_join_1_eigen, y = idx)
    eigen <- Reduce(`+`, eigen)
    
    mds_1_eigen <- mds_1_eigen/length(mds_1_eigen)
    
    eigen <- Reduce(`+`, list(eigen, mds_1_eigen))
    eigen <- eigen/num_partitions
    
    mds_to_return <- list(points = mds_solution, eigen = eigen)
  }
  
  return(mds_to_return)
}
