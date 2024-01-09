# Pruebas con diferentes formas de 
# acelerar cmdscale
#
library(microbenchmark)

n<- 1000
p<- 10
large_lambdas <- c(15,10,5)
main_dim <- length(large_lambdas)
k <-5
Diag_lambda <- diag(c(large_lambdas,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda

D <- as.matrix(dist(X))

# cmdscale
Y0 <- cmdscale(D,k=k, eig=TRUE)
pairs(cbind(X[,1:main_dim],Y0$points))

# cmdscale_svd
source("scripts/pruebas/cmdscale_svd.R")
Ysvd <- cmdscale_svd(D,k=k, eig = TRUE)
pairs(cbind(X[,1:main_dim],Ysvd$points))

pairs(cbind(Y0$points,Ysvd$points))
plot(Y0$eig, Ysvd$eig)
abline(a=0,b=1)

# refund::cmdscale_lanczos
# Pedro: The code is essentially the same as that in 'cmdscale'
# with 'eigen()' replaced by 'mgcv::slanczos()'
# The function mgcv::slanczos calls a C++ function: C_Rlanczos
library(refund)
Y1 <- refund::cmdscale_lanczos(D,k=k)
pairs(cbind(X[,1:main_dim],Y1))

mb1 <- microbenchmark(
  cmdscale(D,k=k),
  cmdscale_svd(D,k=k),
  refund::cmdscale_lanczos(D,k=k),
  times=10)
mb1
plot(mb1)

# It follows that cmdscale_lanczos is much faster! 
# But it uses C code: it calls the mgcv::slanczos function, 
# which calls .C(C_Rlanczos...) which runs in C.
#
# Additionaly, the standard cmdscale is faster than the new cmdscale_svd

library(svd)
source("scripts/pruebas/cmdscale_svd.R")
Y1 <- cmdscale_trlan_svd(D,k=k)
pairs(cbind(X[,1:main_dim],Y1$points))
Y2 <- cmdscale_trlan_eigen(D,k=k)
pairs(cbind(X[,1:main_dim],Y2$points))

mb1 <- microbenchmark(
  #cmdscale(D,k=k),
  cmdscale_trlan_svd(D,k=k),
  cmdscale_trlan_eigen(D,k=k),
  refund::cmdscale_lanczos(D,k=k),
  times=100)
mb1
plot(mb1)



# Rdimtools::do.mds
# From the help: 
# "do.mds performs a classical Multidimensional Scaling (MDS) 
# using Rcpp and RcppArmadillo package to achieve 
# faster performance than cmdscale."
library(Rdimtools)
Y2 <- do.mds(X,ndim=k)
pairs(cbind(X[,1:main_dim],Y2$Y))


# Using irlba::irlba or irlba::partial_eigen
# as in lmds::cmdscale_landmarks
# Pedro: There is an option in irlba() for calling a C++ function: C_IRLB
# This option is used when "fastpath=TRUE"
# 
# Pedro: En el código de "cmdscale_landmarks" sobra la línea 38:
#     dimred <- -t(t(points_inv) %*% ((dist_2lm - rep(mu_n, each = N))/2))
# porque en la línea 39 vuelve a definir dimred como la
# matriz transpuesta de la definida en la línea 38.
# Si se eliminase la línea 38 la función sería un poco más rápida.
# 
# Pedro: La función "cmdscale_landmarks" se puede usar directamente 
# como sustituta de la función cmdscale
library(lmds)
library(irlba)
Y3 <- lmds::cmdscale_landmarks(D,ndim=k,fastpath=TRUE)
pairs(cbind(X[,1:main_dim],Y3))
# redefining cmdscale_landmarks with line 38 commented and replacing
# irlba::partial_eigen by svd::trlan.eigen
source("scripts/pruebas/my_lmds.R") 
Y4 <- my_cmdscale_landmarks(D,ndim=k) # from "my_lmds.R"
pairs(cbind(X[,1:main_dim],Y4))

mb<-  microbenchmark(lmds::cmdscale_landmarks(D,ndim=k,fastpath=TRUE),
                 lmds::cmdscale_landmarks(D,ndim=k),
                 my_cmdscale_landmarks(D,ndim=k),
                 times=10)

mb
plot(mb)

mb <- microbenchmark(
  #cmdscale(D,k=k), 
  refund::cmdscale_lanczos(D,k=k),
  Rdimtools::do.mds(X,ndim=k),
  #lmds::cmdscale_landmarks(D,ndim=k),
  my_cmdscale_landmarks(D,ndim=k),
  cmdscale_trlan_eigen(D,k=k),
  times=100)
mb
plot(mb)
# Conclussions:
# * The fastest one is Rdimtools::do.mds, 
#      then refund::cmdscale_lanczos, 
#      then cmdscale_landmarks (sin línea 38)
# * The quality of lmds::cmdscale_landmarks is very low,
#   where "quality" means ability to recover the 3 main dimensions
# * The quality of Rdimtools::do.mds and refund::cmdscale_lanczos is very high.
#
# So, we will use Rdimtools::do.mds in our simulations.

############################################
#
# Pruebas con diferentes formas de 
# landmark MDS
library(microbenchmark)

n<- 1000000 # 1000000
p<- 100
large_lambdas <- c(15,10,5)
main_dim <- length(large_lambdas)
k <-5
Diag_lambda <- diag(c(large_lambdas,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda

# D <- as.matrix(dist(X))


# lmds::lmds
library(lmds)
nlm <- 1000
Ylm1 <- lmds::lmds(X, ndim=5, num_landmarks = nlm) # no permite fastpath=TRUE
Ismpl <- sample(n,1000)
pairs(cbind(X[Ismpl,1:main_dim],Ylm1[Ismpl,]))

# Rdimtools::do.lmds
library(Rdimtools)
nlm <- 1000
Ylm2 <- Rdimtools::do.lmds(X, ndim=5, npoints = nlm) 
Ylm3 <- Rdimtools::do.lmds(X, ndim=5) 
Ismpl <- sample(n,1000)
pairs(cbind(X[Ismpl,1:main_dim],Ylm2$Y[Ismpl,]))
pairs(cbind(X[Ismpl,1:main_dim],Ylm3$Y[Ismpl,]))

# pivot_mds
source("final_methods/pivot_mds.R")
pivots <- 400
Yp <- pivot_mds(X, pivots=pivots, ndim=5) 
Ismpl <- sample(n,1000)
pairs(cbind(X[Ismpl,1:main_dim],Yp$points[Ismpl,]))


mb.lm <- microbenchmark(
  lmds::lmds(X, ndim=5, num_landmarks = nlm),
  Rdimtools::do.lmds(X, ndim=5, npoints = nlm),
  Rdimtools::do.lmds(X, ndim=5),
  pivot_mds(X, pivots=pivots, ndim=5),
  times=10)
mb.lm 
plot(mb.lm)
# Conclussions:
# * pivot_mds is extreeeeeemly slow!!!!
# * lmds::lmds is extremely slow!!!
# * The quality of lmds::lmss is very low,
#   where "quality" means ability to recover the 3 main dimensions
# * The fastest one is Rdimtools::do.lmds with npoints = 1000, 
#      then Rdimtools::do.lmds with default option for npoints
# * The quality of Rdimtools::do.lmds with default option for npoints is very high.
# * The quality of Rdimtools::do.lmds with npoints=100 is high.

#######################
source("final_methods/divide_conquer_mds.R")
Ydc <- divide_conquer_mds(X, l=400, tie=20, k=5, n_cores=1)
pairs(cbind(X[Ismpl,1:main_dim],Ydc$points[Ismpl,]))

mb.lm <- microbenchmark(
  #lmds::lmds(X, ndim=5, num_landmarks = nlm),
  Rdimtools::do.lmds(X, ndim=5, npoints = nlm),
  Rdimtools::do.lmds(X, ndim=5),
  divide_conquer_mds(X, l=400, tie=20, k=5, n_cores=1),
  lmds::lmds(X, ndim=5, num_landmarks = nlm),
  lmds::lmds(X, ndim=5),
  #pivot_mds(X, pivots=pivots, ndim=5),
  times=10)

plot(mb.lm)
mb.lm
#######################

# Checking results when the columns of X are rotated before doing lmds
rho<-.9
P <- matrix(rho,ncol=dim(X)[2],nrow=dim(X)[2])
diag(P) <- 1
A <- eigen(P)$vectors # this is a rotation matrix
XA <- X %*% A

mb.lm <- microbenchmark(
  #lmds::lmds(X, ndim=5, num_landmarks = nlm),
  Rdimtools::do.lmds(XA, ndim=5, npoints = nlm),
  Rdimtools::do.lmds(XA, ndim=5),
  times=50)

mb.lm
plot(mb.lm)
# conclusion: Rotation of columns of X has no effect on time
YlmA <- Rdimtools::do.lmds(XA, ndim=5)
YlmA <- Rdimtools::do.lmds(XA, ndim=5, npoints = nlm)
Ismpl <- sample(n,1000)
pairs(cbind(XA[Ismpl,1:main_dim],YlmA$Y[Ismpl,]))
pairs(cbind(X[Ismpl,1:main_dim],YlmA$Y[Ismpl,]))
# conclusion: Rotation of columns of X has no effect on quality of the results

#####################################

# Checking Procrustes times using svds
n<- 1000
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

prc <- perform_procrustes(X1, target=X1Ash, X, translation = TRUE)
prcs <- perform_procrustes_svds(X1, target=X1Ash, X, translation = TRUE)
