landmark_mds <- function (x, r, num_landmarks) {
  # Select landmark indexes and compute the distance
  idx_lm <- select_landmarks(x = x, num_landmarks = num_landmarks)

  # Select landmark points
  x_lm <- x[idx_lm, , drop = FALSE]

  # Compute distance with respect those landmark points
  dist_to_lm <- pracma::distmat(X = x_lm, Y = x)

  # Get final configuration
  mds_config <- lmds_config(
    dist_landmark = dist_to_lm,
    idx_landmark = idx_lm,
    x_landmark = x_lm,
    r = r,
    rescale = FALSE
  )
  rownames(mds_config$points) <- rownames(x)
  return(mds_config)
}

lmds_config <- function (dist_landmark, idx_landmark, x_landmark, r, rescale = TRUE){
  # Distance between landmark points and themselves
  dist_lm_lm <- dist_landmark[, idx_landmark, drop = FALSE]^2

  # Obtain number of landmark points
  num_landmarks <- as.integer(nrow(dist_lm_lm))

  # Obtain number of total points
  n_points <- as.integer(ncol(dist_landmark))

  # Double center the matrix
  mu_row_lm <- rowMeans(dist_lm_lm)
  mu_total_lm <- mean(dist_lm_lm)
  x_c <- sweep(dist_lm_lm, 1, mu_row_lm, "-")
  x_dc <- sweep(x_c, 2, mu_row_lm, "-") + mu_total_lm

  # Obtain SVD
  if (nrow(x_dc) > 10) {
    e <- svd::trlan.eigen(-x_dc/2, neig = r)
    ev <- e$d
    evec <- e$u
    if (ncol(evec) != r) {
      e <- eigen(-x_dc/2)
      ev <- e$values[1:r]
      evec <- e$vectors[, 1:r, drop = FALSE]
    }
  } else {
    e <- eigen(-x_dc/2)
    ev <- e$values[1:r]
    evec <- e$vectors[, 1:r, drop = FALSE]
  }

  # Get MDS for all the points
  points_inv <- evec/rep(sqrt(ev), each = num_landmarks)
  mds_config <- (-t(dist_landmark^2 - rep(mu_row_lm, n_points))/2) %*% points_inv
  mds_config <- apply(mds_config, MARGIN = 2, FUN = function(y) y - mean(y))
  mds_config <- mds_config %*% eigen(cov(mds_config))$vectors
  mds <- list()
  mds$points <- mds_config
  mds$eigen <- ev/num_landmarks
  return(mds)
}

select_landmarks <- function (x, num_landmarks) {
  n_row_x <- nrow(x)
  if (num_landmarks > n_row_x) {
    stop("num_landmarks cannot be greater than nrow(x)")
  }

  idx_landmark <- sample.int(n_row_x, num_landmarks)
  return(idx_landmark)
}
