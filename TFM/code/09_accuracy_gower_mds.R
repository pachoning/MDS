source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
source("tools/gower_interpolation_mds.R")


x = matrix(rnorm(5*1000), ncol = 5) %*%diag(c(15,1,1,1,1))
row.names(x) = 1:nrow(x)
s = 5
# Gower MDS
gower_mds = gower.interpolation.mds(
  x = x,
  l = 500,
  s = s,
  metric = "euclidean"
)


# Classical
dist = cluster::daisy(
  x = x,
  metric = "euclidean"
)

mds_classic = stats::cmdscale(
  d = dist, 
  k = s,
  eig = FALSE
)


# Compare
res = compare.methods(
  mds_new_approach = gower_mds$points,
  mds_classical = mds_classic
)

plot(gower_mds$points[,1], res$mds_classical_transformed[,1])
plot(gower_mds$points[,2], res$mds_classical_transformed[,2])
plot(gower_mds$points[,3], res$mds_classical_transformed[,3])
plot(gower_mds$points[,4], res$mds_classical_transformed[,4])
plot(gower_mds$points[,5], res$mds_classical_transformed[,5])







