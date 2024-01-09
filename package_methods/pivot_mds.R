pivot_mds <- function(x, num_pivots, r) {
  n_row_x <- nrow(x)
  if (num_pivots > n_row_x) {
    stop("num_pivots cannot be greater than nrow(x)")
  }

  # Get indexes for pivot points
  idx_pivots <- sample.int(n_row_x, num_pivots)
  x_pivots <- x[idx_pivots, , drop = FALSE]

  # Compute distance matrix
  D <- pracma::distmat(X = x_pivots, Y = x)
  D <- t(D)
  D_sq <- D^2

  # Apply formula
  col_mean <- colMeans(D_sq)
  row_mean <- rowMeans(D_sq)
  Dmat <- -(1/2) * (D_sq - outer(row_mean, col_mean, function(x, y) x + y) + mean(D_sq))

  # Use svd
  if (n_row_x > 10) {
    svd_results <- svd::propack.svd(Dmat, neig = r)
    u_matrix <- svd_results$u
    singular_values <- svd_results$d

    if (ncol(u_matrix) != r) {
      svd_results <- svd(Dmat, nu = r, nv = r)
      u_matrix <- svd_results$u
      singular_values <- svd_results$d
    }
  } else {
    svd_results <- svd(Dmat, nu = r, nv = r)
    u_matrix <- svd_results$u
    singular_values <- svd_results$d
  }

  # Compute MDS
  eigen <- singular_values/sqrt((n_row_x * num_pivots))
  points <- u_matrix %*% diag(sqrt(n_row_x * eigen))
  mds <- list(
    points = points,
    eigen = eigen
  )
  return(mds)
}
