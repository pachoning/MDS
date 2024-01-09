n <- 100000
p <- 100
n1 <- 1000

x = matrix(rnorm(n*p),n,p)

library(pdist)
x.pdist = pdist::pdist(x, indices.A = 1:n1, indices.B = (n1+1):n)
# Converting a pdist object into a traditional distance matrix
D1 <- as.matrix(x.pdist)

library(pracma)
D2 <- pracma::distmat(X=x[1:n1,], Y=x[(n1+1):n,])

max(abs(D1-D2))

library(microbenchmark)
bm <- microbenchmark(
#  pdist::pdist(x, indices.A = 1:n1, indices.B = (n1+1):n),
  as.matrix(pdist::pdist(x, indices.A = 1:n1, indices.B = (n1+1):n)),
  pracma::distmat(X=x[1:n1,], Y=x[(n1+1):n,]),
  dynutils::calculate_distance(x=x[1:n1,], y=x[(n1+1):n,], method="euclidean"),
  times=5
)
bm
plot(bm)
# Concluímos que pracma::distmat es la más rápida

# stats:dist versus pracma::distmat
n1 <- 1000
bm <-  microbenchmark(
  dist(x[1:n1,]),
  pracma::distmat(X=x[1:n1,], Y=x[1:n1,]),
  dynutils::calculate_distance(x=x[1:n1,], y=x[1:n1,], method="euclidean"),
  times=10
)
bm
plot(bm)
# Concluímos que pracma::distmat es la más rápida

