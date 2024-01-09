cmdscale_lanczos <- 
function (d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE) 
{
  if (anyNA(d)) 
    stop("NA values not allowed in 'd'")
  if (is.null(n <- attr(d, "Size"))) {
    if (add) 
      d <- as.matrix(d)
    x <- as.matrix(d^2)
    storage.mode(x) <- "double"
    if ((n <- nrow(x)) != ncol(x)) 
      stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  }
  else {
    rn <- attr(d, "Labels")
    x <- matrix(0, n, n)
    if (add) 
      d0 <- x
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
    if (add) {
      d0[row(x) > col(x)] <- d
      d <- d0 + t(d0)
    }
  }
  n <- as.integer(n)
  if (is.na(n) || n > 46340) 
    stop("invalid value of 'n'")
  if ((k <- as.integer(k)) > n - 1 || k < 1) 
    stop("'k' must be in {1, 2, .. n - 1}")
  x <- scale(t(scale(t(x), scale = FALSE)), scale = FALSE)
  if (add) {
    i2 <- n + (i <- 1L:n)
    Z <- matrix(0, 2L * n, 2L * n)
    Z[cbind(i2, i)] <- -1
    Z[i, i2] <- -x
    Z[i2, i2] <- scale(t(scale(t(2 * d), scale = FALSE)), 
                       scale = FALSE)
    add.c <- max(slanczos(Z, k = 1, kl = 1)$values)
    x <- matrix(double(n * n), n, n)
    non.diag <- row(d) != col(d)
    x[non.diag] <- (d[non.diag] + add.c)^2
    x <- scale(t(scale(t(x), scale = FALSE)), scale = FALSE)
  }
  e <- slanczos(-x/2, k = k)
  ev <- e$values
  evec <- e$vectors
  k1 <- sum(ev > 0)
  if (k1 < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", 
                     k1, k), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }
  points <- evec * rep(sqrt(ev), each = n)
  dimnames(points) <- list(rn, NULL)
  if (eig || x.ret || add) {
    evalus <- e$values
    list(points = points, eig = if (eig) evalus, x = if (x.ret) x, 
         ac = if (add) add.c else 0, GOF = sum(ev)/c(sum(abs(evalus)), 
                                                     sum(pmax(evalus, 0))))
  }
  else points
}
