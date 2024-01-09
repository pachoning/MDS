classical_mds <- function (x, r = 2, return_distance_matrix = FALSE) {
  
  dist_fn = function(x){pracma::distmat(X=x,Y=x)}
  
  # Some sanity checks on the input parameters
  # NA not allowed
  if (anyNA(x)) {
    stop("NA values not allowed in 'x'")
  }
  
  n <- nrow(x)
  
  # Non-numeric, non-integers for r not allowed
  if (is.numeric(r)) {
    if(round(r) != r) {
      stop("'r' be an integer number")
    }
  } 
  
  if (r > n - 1 || r < 1) {
    stop("'r' must be less than nrow(x)}")
  }
  
  rn <- rownames(x)
  
  # Compute distance matrix and double center it
  #dist_matrix <- as.matrix(dist_fn(x))
  dist_matrix <- dist_fn(x)
  dist_2 <- dist_matrix^2
  dist_2 <- scale(t(scale(t(dist_2), scale = FALSE)), scale = FALSE)
  
  # Compute eigen decomposition
  # Avoid this error using svd package. TRLAN is not designed to work with such a small matrix
  if (n >= 10){
    eigen_result <- svd::trlan.eigen(
      -dist_2/2,
      neig = r,
      opts = list(tol = 1e-4)
    )
    eigen_val <- eigen_result$d
    eigen_vec <- eigen_result$u
    # Preventing from getting less columns than r when the eigendecomposition
    # does not converge
    if (ncol(eigen_vec) != r) {
      eigen_result <- eigen(-dist_2/2)
      eigen_val <- eigen_result$values[1:r]
      eigen_vec <- eigen_result$vector[, 1:r, drop = FALSE]
    }
  }else{
    eigen_result <- eigen(-dist_2/2)
    eigen_val <- eigen_result$values[1:r]
    eigen_vec <- eigen_result$vector[, 1:r, drop = FALSE]
  }
  
  # Compute MDS configuration
  points <- eigen_vec * rep(sqrt(eigen_val), each = n)
  rownames(points) <- rn
  mds_result <- list(
    points = points,
    eigen = eigen_val
  )
  
  if (return_distance_matrix) {
    mds_result$distance <- dist_matrix
  }
  
  return(mds_result)
}
