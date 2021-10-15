source("tools/tools_emmanuel.R")

emmanuel_mds <- function(x, l, k, ...){
  
  n <- nrow(x)
  k <- 10
  m0 <- l
  
  ## choose one of the two below:
  # S <- select_m(x, m = m0)
  S <- sample.int(n, size = m0)
  
  d0 <- as.matrix(dist(x[S, ])) # m0 x m0 distance matrix
  MDS <- cmdscale(d0, k = k, eig = TRUE)
  
  q <- k
  
  mds <- MDS$points
  
  ## prepare to use Gower's (1968) formula
  tmds <- t(mds)
  B <- mds %*% tmds
  dB <- diag(B)
  A <- 0.5 * solve(tmds %*% mds) %*% tmds
  
  
  PROJ <- numeric(n * q)
  dim(PROJ) <- c(n, q)
  
  ## the "reduced" (or reference) sample:
  xs <- x[S, ]
  
  iter_loop <- (1:n)[-S] 
  total_iter_loop <- length(iter_loop)
  for (i in iter_loop) {
    tmp <- matrix(x[i, ], m0, ncol(x), byrow = TRUE)
    d <- rowSums((xs - tmp)^2) # assumes Euclidean distance
    ## d <- sqrt(d)
    PROJ[i, ] <- A %*% (dB - d) # avoid squaring d
  }
  PROJ[S, ] <- mds
  
  conf <- list()
  conf$points <- PROJ
  conf$eigen <- rep(1, k)
  conf$GOF <- c(1, 1)
  
  return(conf)

}
