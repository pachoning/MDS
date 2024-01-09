source("final_methods/classical_mds.R")
source("final_methods/landmark_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/interpolation_mds.R")

n <- 8749
p <- 100
large_lambdas <- c(15, 10, 5)
main_dim <- length(large_lambdas)
k <-5
Diag_lambda <- diag(c(large_lambdas^.5,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda
#D <- as.matrix(dist(X))
n_sample <- 100
largest_MDS <- 968

# classical MDS
Y1 <- classical_mds(x = X,k = k)
pairs(cbind(X[,1:main_dim],Y1$points))

# landmark MDS
Y2 <- landmark_mds(x = X,ndim = k,num_landmarks = largest_MDS)
smpl <- sample(dim(X)[1],1000)
pairs(cbind(X[smpl,1:main_dim],Y2$points[smpl,]))

# pivot MDS
Y3 <- pivot_mds(X = X, ndim = k, pivots=largest_MDS)
Y3$eigen
smpl <- sample(dim(X)[1],1000)
pairs(cbind(X[smpl,1:main_dim],Y3$points[smpl,]))

# reduced MDS
Y4 <- reduced_mds(x = X, l = largest_MDS, k = k)
pairs(cbind(X[,1:main_dim],Y4$points))

# divide and conquer MDS
Y5 <- divide_conquer_mds(x = X, l = largest_MDS, tie = n_sample, k = k, n_cores = 1)
pairs(cbind(X[,1:main_dim],Y5$points))

# fast MDS
Y6 <- fast_mds(x = X, l = largest_MDS, s = n_sample, k = k, n_cores = 1)
pairs(cbind(X[,1:main_dim],Y6$points))

# interpolation MDS
Y7 <- interpolation_mds(x = X, l = largest_MDS, k = k, n_cores = 1)
pairs(cbind(X[,1:main_dim],Y7$points))


#################################

n <- 10000
p <- 10
large_lambdas <- c(15, 10, 5)
main_dim <- length(large_lambdas)
k <-5
Diag_lambda <- diag(c(large_lambdas,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda
#D <- as.matrix(dist(X))
n_sample <- 100
largest_MDS <- 1000

library(microbenchmark)
mb <- microbenchmark(
  # classical MDS
  #Y1 <- classical_mds(x = X,k = k),
  # landmark MDS
  Y2 <- landmark_mds(x = X,ndim = k),
  # pivot MDS
  Y3 <- pivot_mds(X = X, ndim = k),
  # reduced MDS
  Y4 <- reduced_mds(x = X, l = largest_MDS, k = k),
  # divide and conquer MDS
  Y5 <- divide_conquer_mds(x = X, l = largest_MDS, tie = n_sample, k = k, n_cores = 1),
  # fast MDS
  Y6 <- fast_mds(x = X, l = largest_MDS, s = n_sample, k = k, n_cores = 1),
  # interpolation MDS
  Y7 <- interpolation_mds(x = X, l = largest_MDS, k = k, n_cores = 1),
  times=10
)
mb
plot(mb)
