# Own version of graphlayouts::layout_with_pmds
# to compute pivot_mds from a data matrix instead of a graph.
#
# X: n\times p data matrix
# pivots: number of pivots
# ndim:	dimensionality of layout (defaults to 2)
pivot_mds <- function(X, pivots = min(dim(X)[1], 100), ndim = 2) {
  n <- dim(X)[1]
  pivs <- sample(1:n, pivots)
  D <- t(pracma::distmat(X=X[pivs,], Y = X))
  cmean <- colMeans(D^2)
  rmean <- rowMeans(D^2)
  Dmat <- -(1/2)*(D^2 - outer(rmean, cmean, function(x, y) x + y) + mean(D^2))
  sl2 <- svd::propack.svd(Dmat, neig = ndim)
  #xy <- Dmat %*% sl2$v[, 1:ndim] %*% diag(sl2$d^(-1/2)) #* sqrt(n/pivots)
  #xy <- sl2$u %*% diag(sl2$d) / (n/pivots)
  eigen <- sl2$d/sqrt((n * pivots))
  xy <- sl2$u %*% diag(sqrt(n * eigen))
  mds <- list(
    points = xy,
    eigen = eigen
  )
  return(mds)
}
