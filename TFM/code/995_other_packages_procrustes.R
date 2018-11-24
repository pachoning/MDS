data(varespec)
vare.dist <- vegdist(wisconsin(varespec))
library(MASS)  ## isoMDS
mds.null <- isoMDS(vare.dist, tol=1e-7)
mds.alt <- isoMDS(vare.dist, initMDS(vare.dist), maxit=200, tol=1e-7)
vare.proc <- vegan::procrustes(mds.alt$points, mds.null$points, symmetric = TRUE)
vare.proc$scale
sss = pracma::procrustes(mds.alt$points, mds.null$points)
head(sss$P)
head(mds.null$points %*% vare.proc$rotation)

head(mds.null$points %*% vare.proc$rotation - mds.alt$points)
head(vare.proc$Yrot- mds.alt$points)
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)

?residuals

rotate_matrix_divide = procrustes(
  mds_x,
  mds_algorithm$mds
)