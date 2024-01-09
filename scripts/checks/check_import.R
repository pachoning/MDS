library(Rcpp)
library(microbenchmark)
library(corpcor)
sourceCpp("scripts/pruebas/prueba.cpp")

n<- 100000
n1 <- 100
p<- 10
large_lambdas <- c(15,10,5)
main_dim <- length(large_lambdas)
k <-5
Diag_lambda <- diag(c(large_lambdas,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda

rho<-.9
P <- matrix(rho,ncol=dim(X)[2],nrow=dim(X)[2])
diag(P) <- 1
A <- eigen(P)$vectors # this is a rotation matrix


X1 <- X[1:n1,]
shift <- rnorm(dim(X1)[2],m=5)
X1Ash <- X1 %*% A + matrix(shift,nrow=dim(X1)[1],ncol=dim(X1)[2],byrow=TRUE)
XAsh <- X %*% A + matrix(shift,nrow=dim(X)[1],ncol=dim(X)[2],byrow=TRUE)


result <- getEigenValues(X)

pairs(cbind(X[, 1:3], result$leftSingularVectors[, 1:3]))
svd_result <- svd(X)
result$singularValues
svd_result$d


max(abs(result$leftSingularVectors) - abs(svd_result$u))
dim(result$leftSingularVectors)
dim(svd_result$u)


max(abs(result$rightSingularVectors) - abs(svd_result$v))
min(abs(result$rightSingularVectors) - abs(svd_result$v))

mb <- microbenchmark(
  fast.svd(XAsh),
  svd(XAsh),
  times=50
)
