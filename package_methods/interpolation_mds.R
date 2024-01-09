source("package_methods/classical_mds.R")
source("tools/procrustes.R")

get_partitions_for_interpolation <- function(n, n_obs, l, r) {

  if (l <= r) {
    stop("l must be greater than r")
  }

  if (n <= l) {
    p <- 1
  } else {
    p <- 1 + ceiling((n - l)/n_obs)
    n_last <- n - (l + (p - 2) * n_obs)
  }

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (p > 1) {
    first_part <- permutation[1:l]
    middle_part <- permutation[(l + 1):(n - n_last)]
    last_part <- permutation[(n - n_last + 1):n]

    list_index <- split(middle_part, 1:(p - 2))
    names(list_index) <- NULL
    list_index[[p - 1]] <- list_index[[1]]
    list_index[[1]] <- first_part
    list_index[[p]] <- last_part

  } else {
    list_index <- list(permutation)
  }

  return(list_index)
}

get_P_matrix <- function(n_row) {

  identity_matrix <- diag(x = 1, nrow = n_row, ncol = n_row)
  one_vector <- matrix(data = 1, nrow = n_row, ncol = 1)
  P <- identity_matrix - 1/n_row * one_vector %*% t(one_vector)
  return(P)

}

interpolation_mds_main <- function(idx, x, data_1, x_1, n_row_1, q_vector, x_1__s_1__inv) {

  # Filter the matrix
  x_other <- x[idx, , drop = FALSE]
  n_row_other <- nrow(x_other)
  n_row_1 <- nrow(data_1)

  # Get A matrix
  A <- pracma::distmat(X = x_other, Y = data_1)

  # Get delta matrix
  A_sq <- A^2

  # One vecotr
  one_vector_other <- matrix(data = 1, nrow = n_row_other, ncol = 1)

  # Get coordinates
  x_2 <- 1/(2 * n_row_1) * (one_vector_other %*% t(q_vector) - A_sq) %*% x_1__s_1__inv

  return(x_2)
}

interpolation_mds <- function(x, l, r, n_cores) {

  n <- nrow(x)
  n_row_partition <- l
  indexes_partition <- get_partitions_for_interpolation(n = n, n_obs = n_row_partition, l = l, r = r)
  num_partitions <- length(indexes_partition)

  if (num_partitions <= 1) {
    # It is possible to run MDS directly
    mds <- classical_mds(x = x, r = r)
    points <- mds$points
    eigen_v <- mds$eigen/n
    list_to_return <- list(points = points, eigen = eigen_v)
  } else {
    # Get the first group 
    n_row_1 <- length(indexes_partition[[1]])

    # Obtain MDS for the first group
    data_1 <- x[indexes_partition[[1]], , drop = FALSE]
    mds_eig <- classical_mds(x = data_1, r = r, return_distance_matrix = TRUE)
    distance_matrix <- as.matrix(mds_eig$distance)
    X_1 <- mds_eig$points
    eigen_v <- mds_eig$eigen/nrow(X_1)

    # Get P matrix
    P <- get_P_matrix(n_row = n_row_1)
    Q <- -1/2 * P %*% distance_matrix^2 %*% t(P)
    q_vector <- diag(Q)
    S <- 1 / (n_row_1-1) * t(X_1) %*% X_1
    x_1__s_1__inv <- X_1 %*% solve(S)

    # Calculations needed to do Gower interpolation
    mds_others <- parallel::mclapply(
      indexes_partition[2:num_partitions],
      interpolation_mds_main,
      x = x,
      data_1 = data_1,
      x_1 = X_1,
      n_row_1 = n_row_1,
      q_vector = q_vector,
      x_1__s_1__inv = x_1__s_1__inv,
      mc.cores = n_cores
    )

    mds_points <- matrix(data = NA, nrow = n, ncol = r)
    mds_points[1:n_row_1, ] <- X_1
    mds_points[(n_row_1 + 1):n, ] <- do.call(rbind, mds_others)
    idx_all <- do.call(c, indexes_partition)
    idx_all <- order(idx_all)
    mds_points <- mds_points[idx_all, , drop = FALSE]
    mds_points <- apply(mds_points, MARGIN = 2, FUN = function(y) y - mean(y))
    mds_points <- mds_points %*% eigen(cov(mds_points))$vectors
    list_to_return <- list(points = mds_points, eigen = eigen_v)
  }

  return(list_to_return)
}
