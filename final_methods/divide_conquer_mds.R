source("final_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer <- function(n, l, tie, k) {
  
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


filter_matrix <- function(matrix, idx) {
  return(matrix[idx, , drop = FALSE])
}


join_matrices <- function(m1, m2) {
  return(rbind(m1, m2))
}

split_matrix <- function(matrix, num_points) {
  x_1 <- matrix[1:num_points, , drop = FALSE]
  x_2 <- matrix[(num_points + 1):nrow(matrix), , drop = FALSE]
  return(list(x_1, x_2))
}

#main_divide_conquer_mds <- function(idx, x, x_sample_1, k, original_mds_sample_1, dist_fn, ...) {
main_divide_conquer_mds <- function(idx, x, x_sample_1, k, original_mds_sample_1, dist_fn) {  
  # Filter the matrix
  x_filtered <- filter_matrix(matrix = x, idx = idx)
  
  # Join with the sample from the first partition
  x_join_sample_1 <- join_matrices(m1 = x_sample_1, m2 = x_filtered)
  
  # Perform MDS
  #mds_all <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
  mds_all <- classical_mds(x = x_join_sample_1, k = k, dist_fn = dist_fn)
  mds_points <- mds_all$points
  mds_eigen <- mds_all$eigen
  mds_GOF <- mds_all$GOF
  
  # Split MDS into two parts: the first part is the MDS for the sampling points on the
  # first partition and the second part is the MDS for the points of the partition
  mds_split <- split_matrix(matrix = mds_points, num_points = nrow(x_sample_1))
  mds_sample_1 <- mds_split[[1]]
  mds_partition <- mds_split[[2]]
  mds_procrustes <- perform_procrustes(x = mds_sample_1, 
                                       target = original_mds_sample_1, 
                                       matrix_to_transform = mds_partition, 
                                       translation = FALSE)
  
  mds_eigen <- mds_eigen/length(idx)
  mds_GOF <- mds_GOF*length(idx)
  
  return(list(points = mds_procrustes, eigen = mds_eigen, GOF = mds_GOF))
}

divide_conquer_mds <- function(x, l, tie, k, dist_fn = stats::dist, n_cores, ...) {
  
  n_row_x <- nrow(x)
  
  if (n_row_x <= l) {
    #mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    mds_to_return$eigen <- mds_to_return$eigen/n_row_x
    mds_to_return$GOF <- mds_to_return$GOF
    
  } else {
    
    mds_matrix <- matrix(data = NA, nrow = n_row_x, ncol = k)
    
    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, tie = tie, k = k)
    num_partitions <- length(idx)
    length_1 <- length(idx[[1]])
    
    # Perform MDS for the first partition
    x_1 <- x[idx[[1]], , drop = FALSE]
    #mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, ...)
    mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/length_1
    mds_1_GOF <- mds_1$GOF*length_1
    
    # Take a sample from the first partition
    sample_1 <- sample(x = length_1, size = tie, replace = FALSE)
    x_sample_1 <- x_1[sample_1, ,drop = FALSE]
    mds_sample_1 <- mds_1_points[sample_1, , drop = FALSE]
    
    #mds_others_results <- parallel::mclapply(idx[2:num_partitions],
    #                      main_divide_conquer_mds,
    #                      x = x, 
    #                      x_sample_1 = x_sample_1, 
    #                      k = k, 
    #                      original_mds_sample_1 = mds_sample_1, 
    #                      dist_fn = dist_fn,
    #                     mc.cores = n_cores,
    #                      ...)
    
    mds_others_results <- parallel::mclapply(idx[2:num_partitions],
                                             main_divide_conquer_mds,
                                             x = x,
                                             x_sample_1 = x_sample_1,
                                             k = k,
                                             original_mds_sample_1 = mds_sample_1,
                                             dist_fn = dist_fn,
                                             mc.cores = n_cores)
    
    
    # Obtain points
    mds_others_points <- do.call(rbind, parallel::mclapply(mds_others_results, function(x) x$points, mc.cores = n_cores))
    mds_matrix[1:length_1, ] <- mds_1_points
    mds_matrix[(length_1 + 1):n_row_x, ] <- mds_others_points
    order_idx <- do.call(c, idx)
    order_idx <- order(order_idx)
    mds_matrix <- mds_matrix[order_idx, , drop = FALSE]
    mds_matrix <- apply(mds_matrix, MARGIN = 2, FUN = function(x) x - mean(x))
    mds_matrix <- mds_matrix %*% eigen(cov(mds_matrix))$vectors
    
    # Obtain eigenvalues
    eigen <- parallel::mclapply(mds_others_results, function(x) x$eigen, mc.cores = n_cores)
    eigen[[num_partitions]] <- mds_1_eigen
    eigen <- Reduce(`+`, eigen)
    eigen <- eigen/num_partitions
    
    # Obtain GOF
    GOF <- parallel::mclapply(mds_others_results, function(x) x$GOF, mc.cores = n_cores)
    GOF[[num_partitions]] <- mds_1_GOF
    GOF <- Reduce(`+`, GOF)
    GOF <- GOF/n_row_x
    mds_to_return <- list(points = mds_matrix, eigen = eigen, GOF = GOF)
  }
  
  return(mds_to_return)
}
