source("package_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_divide_conquer <- function(n, l, c_points, r) {

  if (l - c_points <= 0) {
    stop("l must be greater than c_points")
  } else if(l - c_points < c_points) {
    stop("l-c_points must be greater than c_points")
  } else if(l - c_points <= r){
    stop("l-c_points must be greater than r")
  }

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (n<=l) {
    list_indexes <- list(permutation)
  } else {
    min_sample_size <- max(r + 2, c_points)
    p <- 1 + ceiling((n - l)/(l - c_points))
    last_partition_sample_size <- n - (l + (p-2) * (l - c_points))

    if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
      p <- p - 1
      last_partition_sample_size <- n - (l + (p-2) * (l-c_points))
    }

    first_parition <- permutation[1:l]
    last_partition <- permutation[(n - last_partition_sample_size + 1):n]
    list_indexes <- split(x = permutation[(l + 1):(n - last_partition_sample_size)], f = 1:(p - 2))
    names(list_indexes) <- NULL
    list_indexes[[p - 1]] <- list_indexes[[1]]
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

main_divide_conquer_mds <- function(idx, x, x_sample_1, r, original_mds_sample_1) {  
  # Filter the matrix
  x_filtered <- filter_matrix(matrix = x, idx = idx)

  # Join with the sample from the first partition
  x_join_sample_1 <- join_matrices(m1 = x_sample_1, m2 = x_filtered)

  # Perform MDS
  mds_all <- classical_mds(x = x_join_sample_1, r = r)
  mds_points <- mds_all$points
  mds_eigen <- mds_all$eigen

  # Split MDS into two parts: the first part is the MDS for the sampling points on the
  # first partition and the second part is the MDS for the points of the partition
  mds_split <- split_matrix(matrix = mds_points, num_points = nrow(x_sample_1))
  mds_sample_1 <- mds_split[[1]]
  mds_partition <- mds_split[[2]]
  mds_procrustes <- perform_procrustes(
    x = mds_sample_1,
    target = original_mds_sample_1,
    matrix_to_transform = mds_partition,
    translation = FALSE
  )

  mds_eigen <- mds_eigen/length(idx)

  return(list(points = mds_procrustes, eigen = mds_eigen))
}

divide_conquer_mds <- function(x, l, c_points, r, n_cores) {
  n_row_x <- nrow(x)
  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, r = r)
    mds_to_return$eigen <- mds_to_return$eigen/n_row_x
  } else {
    mds_matrix <- matrix(data = NA, nrow = n_row_x, ncol = r)

    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, c_points = c_points, r = r)
    num_partitions <- length(idx)
    length_1 <- length(idx[[1]])

    # Perform MDS for the first partition
    x_1 <- x[idx[[1]], , drop = FALSE]
    mds_1 <- classical_mds(x = x_1, r = r)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/length_1

    # Take a sample from the first partition
    sample_1 <- sample(x = length_1, size = c_points, replace = FALSE)
    x_sample_1 <- x_1[sample_1, ,drop = FALSE]
    mds_sample_1 <- mds_1_points[sample_1, , drop = FALSE]

    mds_others_results <- parallel::mclapply(
      idx[2:num_partitions],
      main_divide_conquer_mds,
      x = x,
      x_sample_1 = x_sample_1,
      r = r,
      original_mds_sample_1 = mds_sample_1,
      mc.cores = n_cores
    )

    # Obtain points
    mds_others_points <- do.call(rbind, parallel::mclapply(mds_others_results, function(x) x$points, mc.cores = n_cores))
    mds_matrix[1:length_1, ] <- mds_1_points
    mds_matrix[(length_1 + 1):n_row_x, ] <- mds_others_points
    order_idx <- do.call(c, idx)
    order_idx <- order(order_idx)
    mds_matrix <- mds_matrix[order_idx, , drop = FALSE]
    mds_matrix <- apply(mds_matrix, MARGIN = 2, FUN = function(x) x - mean(x))
    cov_mds_matrix <- cov(mds_matrix)
    svd_mds_matrix <- eigen(cov_mds_matrix)
    eigenvectors_cov <- svd_mds_matrix$vectors
    mds_matrix <- mds_matrix %*% eigenvectors_cov

    # Obtain eigenvalues
    eigen <- parallel::mclapply(mds_others_results, function(x) x$eigen, mc.cores = n_cores)
    eigen[[num_partitions]] <- mds_1_eigen
    eigen <- Reduce(`+`, eigen)
    eigen <- eigen/num_partitions

    mds_to_return <- list(points = mds_matrix, eigen = eigen)
  }

  return(mds_to_return)
}
