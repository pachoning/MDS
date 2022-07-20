library(lmds)
my_lmds <- function(x, k, l,...) {
  result <- lmds::lmds(x = x, ndim = k, num_landmarks = l)
  conf <- list()
  conf$points <- result
  conf$eigen <- rep(1, k)
  conf$GOF <- c(1, 1)
  return(conf)
}