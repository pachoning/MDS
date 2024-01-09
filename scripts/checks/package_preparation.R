source("package_methods/classical_mds.R")
source("package_methods/landmark_mds.R")
source("package_methods/interpolation_mds.R")
source("package_methods/reduced_mds.R")
source("package_methods/pivot_mds.R")
source("package_methods/divide_conquer_mds.R")
source("package_methods/fast_mds.R")

n <- 8749
p <- 100
large_lambdas <- c(15, 10, 5)
main_dim <- length(large_lambdas)
r <-5
Diag_lambda <- diag(c(large_lambdas^.5,rep(1,(p-main_dim))))
X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Diag_lambda
n_sample <- 100
largest_MDS <- 968

# classical MDS
Y1 <- classical_mds(x = X,r = r)
pairs(cbind(X[,1:main_dim],Y1$points))

# landmark MDS
Y2 <- landmark_mds(x = X, r = r, num_landmarks = largest_MDS)
smpl <- sample(dim(X)[1], 1000)
pairs(cbind(X[smpl,1:main_dim], Y2$points[smpl,]))
Y2$eigen

# interpolation MDS
Y3 <- interpolation_mds(x = X, l = largest_MDS, r = r, n_cores = 1)
pairs(cbind(X[smpl,1:main_dim], Y3$points[smpl,]))
Y3$eigen

# reduced MDS
Y4 <- reduced_mds(x = X, l = largest_MDS, r = r, n_cores= 1)
pairs(cbind(X[smpl,1:main_dim], Y4$points[smpl,]))
Y4$eigen

# pivot MDS
Y5 <- pivot_mds(x = X, num_pivots = largest_MDS, r = r)
pairs(cbind(X[smpl,1:main_dim], Y5$points[smpl,]))
Y5$eigen

# divide-and-conquer MDS
Y6 <- divide_conquer_mds(x = X, l = largest_MDS, c_points = 5 * r, r = r, n_cores = 1)
pairs(cbind(X[smpl,1:main_dim], Y6$points[smpl,]))
Y6$eigen

# fast MDS
Y7 <- fast_mds(x = X, l = largest_MDS, s_points  = 5 * r, r = r, n_cores = 1)
pairs(cbind(X[smpl,1:main_dim], Y7$points[smpl,]))
Y7$eigen
