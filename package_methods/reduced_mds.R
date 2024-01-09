source("package_methods/classical_mds.R")
reduced_mds <- function(x, l, r, n_cores) {
  n <- nrow(x)

  # Take a random sample
  idx <- sample.int(n, size = l)

  # Get MDS for the random sample
  x_initial <- x[idx, , drop = FALSE]
  x_initial_mds <- classical_mds(x = x_initial, r = r)
  mds_initial <- x_initial_mds$points

  # Prepare to use Gower's formula
  t_mds_initial <- t(mds_initial)
  Q <- mds_initial %*% t_mds_initial
  q <- diag(Q)
  A <- 0.5 * solve(t_mds_initial %*% mds_initial) %*% t_mds_initial

  mds_config <- numeric(n * r)
  dim(mds_config) <- c(n, r)

  # Get config for the remaining points
  idx_remaining <- (1:n)[-idx]
  mds_remaining_list <- parallel::mclapply(
    idx_remaining,
    partial_reduced_mds,
    x = x,
    l = l,
    x_initial = x_initial,
    A = A,
    q = q,
    mc.cores = n_cores
  )
  mds_remaining <- do.call(rbind, mds_remaining_list)
  # Populate final mds with initial mds
  mds_config[idx, ] <- mds_initial
  # Populate final mds with remaining mds
  mds_config[idx_remaining, ] <- mds_remaining
  # Normalise
  mds_config <- apply(mds_config, MARGIN = 2, FUN = function(y) y - mean(y))
  mds_config <- mds_config %*% base::eigen(stats::cov(mds_config))$vectors
  # Return final object
  conf <- list()
  conf$points <- mds_config
  conf$eigen <- x_initial_mds$eigen/l
  return(conf)
}

partial_reduced_mds <- function(i, x, l, x_initial, A, q) {
  current_point <-  matrix(x[i, , drop = FALSE], l, ncol(x), byrow = TRUE)
  dist_points <- rowSums((x_initial - current_point)^2)
  result <- t(A %*% (q - dist_points))
  return(result)
}
