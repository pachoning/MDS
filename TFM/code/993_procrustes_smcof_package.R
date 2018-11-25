# https://rdrr.io/cran/smacof/man/Procrustes.html
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")

## artificial example:
X <- matrix(c(1, -1, -1, 1, 2, 2, -2, -2), ncol = 2)
Y <- matrix(c(0.07, 0.93, 1.93, 1.07, 2.62, 3.12, 1.38, 0.88), ncol = 2)
op <- par(mfrow = c(1,2))
plot(X[,1], X[,2], xlim = c(-3, 3), ylim = c(-2, 3.5), asp = 1, xlab = "", ylab = "")
rect(-1, -2, 1, 2)
points(Y[,1], Y[,2], xlim = c(-3, 3), col = "gray")
polygon(Y[,1], Y[,2], border = "gray")
fitp <- Procrustes(X, Y)
plot(fitp$Yhat[,1], fitp$Yhat[,2], col = "red", xlim = c(-3, 3), ylim = c(-2, 3.5), 
     asp = 1, xlab = "", ylab = "")
polygon(fitp$Yhat[,1], fitp$Yhat[,2], border = "red")
par(op)

## MDS example:
eastD <- sim2diss(EW_eng$east)
attr(eastD, "Labels") <- abbreviate(attr(eastD, "Labels"))
fit.east <- mds(eastD, type = "ordinal")
westD <- sim2diss(EW_eng$west)
attr(westD, "Labels") <- abbreviate(attr(westD, "Labels"))
fit.west <- mds(westD, type = "ordinal", init = torgerson(eastD))

fit.proc <- Procrustes(fit.east$conf, fit.west$conf)
fit.proc$Yhat

## Checking the outputs
head(fit.east$conf)
head(fit.proc$dilation * as.matrix(fit.west$conf) %*% as.matrix(fit.proc$rotation) + fit.proc$translation)
head(fit.proc$Yhat)


## Configuration plots; Procrustes plots.
plot(fit.east, main = "MDS East Germany")   ## MDS plot East Germany
plot(fit.west, main = "MDS West Germany")   ## MDS plot West Germany

## Procrustes configurations (X and Yhat)
plot(fit.proc, ylim = c(-1, 1),  col.X = "cadetblue", col.Yhat = "brown", pch = 19, 
     legend = list(pos = "topleft", labels = c("East Germany", "West Germany"))) 

## Procrustes transformations (Y and Yhat)
plot(fit.proc, plot.type = "transplot", length = 0.05, ylim = c(-1,1), 
     legend = list(pos = "bottomright", 
                   labels = c("West Germany (untransformed)", "West Germany (transformed)")))


## Very random example
set.seed(12345)
n_obs = 10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)


dist_matrix = daisy(
  x = x,
  metric = "euclidean"
)

mds_classical = stats::cmdscale(
  d = dist_matrix,
  k = 2
)
  

groups = sample(x = 2, size = nrow(x), replace = TRUE)

mds_divide_conquer_mds = divide_conquer_mds(
  x = x,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)

head(mds_divide_conquer_mds$ls_positions[[1]])
head(x)
head(mds_divide_conquer_mds$mds)
head(mds_classical)

procrustes_result2 =  smacof::Procrustes(
  X = mds_classical[1:20,], 
  Y = mds_divide_conquer_mds$mds[1:20,]
)

head(mds_divide_conquer_mds$mds %*% procrustes_result2$rotation)
head(mds_classical)

plot(
  procrustes_result2, 
  ylim = c(-3, 3),  
  col.X = "cadetblue", 
  col.Yhat = "brown", 
  pch = 2, 
  legend = list(
    pos = "topleft"
  )
) 



x_perm = x[nrow(x):1, ]


dist_matrix = daisy(
  x = x_perm,
  metric = "euclidean"
)


mds_classical2 = stats::cmdscale(
  d = dist_matrix,
  k = 2
)

mds_classical2 = mds_classical2[nrow(x):1, ]

head(mds_classical)
head(mds_classical2)

proc2 = Procrustes(mds_classical, mds_classical2)
head(mds_classical)
head(proc2$Yhat)
