classical_mds <- function (x, k = 2, return_distance_matrix = FALSE) {
  
  dist_fn = function(x){pracma::distmat(X=x,Y=x)}
  
  # Some sanity checks on the input parameters
  # NA not allowed
  if (anyNA(x)) {
    stop("NA values not allowed in 'x'")
  }
  
  n <- nrow(x)
  
  # Non-numeric, non-integers for k not allowed
  if (is.numeric(k)) {
    if(round(k) != k) {
      stop("'k' be an integer number")
    }
  } 
  
  if (k > n - 1 || k < 1) {
    stop("'k' must be less than nrow(x)}")
  }
  
  rn <- rownames(x)
  
  # Compute distance matrix and double center it
  #dist_matrix <- as.matrix(dist_fn(x))
  dist_matrix <- dist_fn(x)
  dist_2 <- dist_matrix^2
  dist_2 <- scale(t(scale(t(dist_2), scale = FALSE)), scale = FALSE)
  
  # Compute eigen decomposition
  
  # Avoid this error using svd package
  # TRLAN is not designed to work with such a small matrix
  if (n >= 10){
    eigen_result <- svd::trlan.eigen(
      -dist_2/2,
      neig = k,
      opts = list(tol = 1e-4)
    )
    eigen_val <- eigen_result$d
    eigen_vec <- eigen_result$u
    # Preventing from getting less columns than k when the eigendecomposition
    # does not converge
    if (ncol(eigen_vec) != k) {
      eigen_result <- eigen(-dist_2/2)
      eigen_val <- eigen_result$values[1:k]
      eigen_vec <- eigen_result$vector[, 1:k, drop = FALSE]
    }
  }else{
    eigen_result <- eigen(-dist_2/2)
    eigen_val <- eigen_result$values[1:k]
    eigen_vec <- eigen_result$vector[, 1:k, drop = FALSE]
  }
  
  #eigen_result <- eigen(-dist_2/2)
  #eigen_val <- eigen_result$values[1:k]
  #eigen_vec <- eigen_result$vector[, 1:k]
  
  # Take positive eigenvalues
  #k1 <- sum(eigen_val > 0)
  #if (k1 < k) {
  #  warning(gettextf("only %d of the first %d eigenvalues are > 0", k1, k), domain = NA)
  #  eigen_vec <- eigen_vec[, eigen_val > 0, drop = FALSE]
  #  eigen_val <- eigen_val[eigen_val > 0]
  #}
  
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
