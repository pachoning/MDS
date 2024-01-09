# Version of stats::cmdscale using svd instead of eigen, 
# because we usually want just the first k dimensions.
# eigen compute all eigenvalues and eigenvectors.
# svd can compute only the first k of them.
cmdscale_svd <- 
  function (d, k = 2, eig = FALSE, x.ret = FALSE) 
  {
    if (anyNA(d)) 
      stop("NA values not allowed in 'd'")
    if (is.null(n <- attr(d, "Size"))) {
      x <- as.matrix(d^2)
      storage.mode(x) <- "double"
      if ((n <- nrow(x)) != ncol(x)) 
        stop("distances must be result of 'dist' or a square matrix")
      rn <- rownames(x)
    }
    else {
      rn <- attr(d, "Labels")
      x <- matrix(0, n, n)
      x[row(x) > col(x)] <- d^2
      x <- x + t(x)
    }
    n <- as.integer(n)
    if (is.na(n) || n > 46340) 
      stop(gettextf("invalid value of %s", "'n'"), domain = NA)
    if ((k <- as.integer(k)) > n - 1 || k < 1) 
      stop("'k' must be in {1, 2, ..  n - 1}")
    x <- scale(t(scale(t(x), scale = FALSE)), scale = FALSE)
    #e <- eigen(-x/2, symmetric = TRUE)
    e <- svd(-x/2, nu=k, nv=k)
    #ev <- e$values[seq_len(k)]
    ev <- e$d[1:k]
    #evec <- e$vectors[, seq_len(k), drop = FALSE]
    evec <- e$u
    k1 <- sum(ev > 0)
    if (k1 < k) {
      warning(gettextf("only %d of the first %d eigenvalues are > 0", 
                       k1, k), domain = NA)
      evec <- evec[, ev > 0, drop = FALSE]
      ev <- ev[ev > 0]
    }
    points <- evec * rep(sqrt(ev), each = n)
    dimnames(points) <- list(rn, NULL)
    evalus <- e$d
    return(list(points = points, eig = if (eig) evalus, x = if (x.ret) x, 
                GOF = sum(ev)/c(sum(abs(evalus)), sum(pmax(evalus, 0)))))
  }

cmdscale_trlan_svd <- 
  function (d, k = 2, eig = FALSE, x.ret = FALSE) 
  {
    require(svd)
    if (anyNA(d)) 
      stop("NA values not allowed in 'd'")
    if (is.null(n <- attr(d, "Size"))) {
      x <- as.matrix(d^2)
      storage.mode(x) <- "double"
      if ((n <- nrow(x)) != ncol(x)) 
        stop("distances must be result of 'dist' or a square matrix")
      rn <- rownames(x)
    }
    else {
      rn <- attr(d, "Labels")
      x <- matrix(0, n, n)
      x[row(x) > col(x)] <- d^2
      x <- x + t(x)
    }
    n <- as.integer(n)
    if (is.na(n) || n > 46340) 
      stop(gettextf("invalid value of %s", "'n'"), domain = NA)
    if ((k <- as.integer(k)) > n - 1 || k < 1) 
      stop("'k' must be in {1, 2, ..  n - 1}")
    x <- scale(t(scale(t(x), scale = FALSE)), scale = FALSE)
    #e <- eigen(-x/2, symmetric = TRUE)
    e <- trlan.svd(-x/2, neig=k)
    #ev <- e$values[seq_len(k)]
    ev <- e$d[1:k]
    #evec <- e$vectors[, seq_len(k), drop = FALSE]
    evec <- e$u
    k1 <- sum(ev > 0)
    if (k1 < k) {
      warning(gettextf("only %d of the first %d eigenvalues are > 0", 
                       k1, k), domain = NA)
      evec <- evec[, ev > 0, drop = FALSE]
      ev <- ev[ev > 0]
    }
    points <- evec * rep(sqrt(ev), each = n)
    dimnames(points) <- list(rn, NULL)
    evalus <- e$d
    return(list(points = points, eig = if (eig) evalus, x = if (x.ret) x, 
                GOF = sum(ev)/c(sum(abs(evalus)), sum(pmax(evalus, 0)))))
  }

cmdscale_trlan_eigen <- 
  function (d, k = 2, eig = FALSE, x.ret = FALSE) 
  {
    require(svd)
    if (anyNA(d)) 
      stop("NA values not allowed in 'd'")
    if (is.null(n <- attr(d, "Size"))) {
      x <- as.matrix(d^2)
      storage.mode(x) <- "double"
      if ((n <- nrow(x)) != ncol(x)) 
        stop("distances must be result of 'dist' or a square matrix")
      rn <- rownames(x)
    }
    else {
      rn <- attr(d, "Labels")
      x <- matrix(0, n, n)
      x[row(x) > col(x)] <- d^2
      x <- x + t(x)
    }
    n <- as.integer(n)
    if (is.na(n) || n > 46340) 
      stop(gettextf("invalid value of %s", "'n'"), domain = NA)
    if ((k <- as.integer(k)) > n - 1 || k < 1) 
      stop("'k' must be in {1, 2, ..  n - 1}")
    x <- scale(t(scale(t(x), scale = FALSE)), scale = FALSE)
    #e <- eigen(-x/2, symmetric = TRUE)
    e <- trlan.eigen(-x/2, neig=k)
    #ev <- e$values[seq_len(k)]
    ev <- e$d[1:k]
    #evec <- e$vectors[, seq_len(k), drop = FALSE]
    evec <- e$u
    k1 <- sum(ev > 0)
    if (k1 < k) {
      warning(gettextf("only %d of the first %d eigenvalues are > 0", 
                       k1, k), domain = NA)
      evec <- evec[, ev > 0, drop = FALSE]
      ev <- ev[ev > 0]
    }
    points <- evec * rep(sqrt(ev), each = n)
    dimnames(points) <- list(rn, NULL)
    evalus <- e$d
    return(list(points = points, eig = if (eig) evalus, x = if (x.ret) x, 
                GOF = sum(ev)/c(sum(abs(evalus)), sum(pmax(evalus, 0)))))
  }

