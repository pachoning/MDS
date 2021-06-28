classical_mds <- function(x, k, dist_fn, return_distance_matrix = FALSE, ...) {

  if (!is.function(dist_fn)) {
    stop("dist_fn must be a function")
  } else if(any(is.na(x))) {
    stop("There are some NA values in the data. Please remove them")
  }

  mds <- list()
  dist_matrix <- dist_fn(x, ...)
  mds_result <- stats::cmdscale(d = dist_matrix, k = k, eig = TRUE)

  mds$points <- mds_result$points
  mds$eigen <- mds_result$eig[1:k]
  mds$GOF <- mds_result$GOF

  if (return_distance_matrix) {
    mds$distance <- as.matrix(dist_matrix)
  }

  return(mds)
}
