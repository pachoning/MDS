## function to get the diagnostic to assess the convergence of the reduced MDS
## d: distance matrix (must be square, symmetric, with n rows and n columns)
## n: number of rows and of columns of d
getDiag <- function(d, n)
{
  ##d <- .Call(stats:::C_DoubleCentre, d^2)
  ##d <- -d/2
  e <- eigen(d, symmetric = TRUE)
  vars <- e$values
  vars[vars < 0] <- 0
  z <- e$vectors * rep(sqrt(vars), each = n)
  dimsel <- which(vars > 0)
  npairs <- (n * (n - 1) / 2)
  dd <- as.dist(d) - dist(z[, dimsel, drop = FALSE])
  sd <- sqrt(sum(dd^2) / npairs)
  m <- sum(dd) / npairs
  sd / m # CV
}

## function to find the "optimal" value of m
opt.m <- function(d0, S, quiet = FALSE)
{
  m <- 3L
  i <- 1:m
  crit <- numeric(m0 <- length(S))
  d0 <- as.matrix(d0)
  crit[3L] <- getDiag(d0[i, i], 3L)
  while (abs(crit[m] - crit[m - 1L]) > 1e-6 && m < m0) {
    if (!quiet) cat("\rm =", m)
    m <- m + 1L
    i <- 1:m
    crit[m] <- getDiag(d0[i, i], m)
  }
  if (!quiet) cat("\n")
  crit[1:m]
}