source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")


set.seed(12345)

# 01 Very random example ----
n_obs = 10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)

# Classical MDS
dist_matrix = daisy(
  x = x,
  metric = "euclidean"
)


mds_classical = stats::cmdscale(
  d = dist_matrix,
  k = 2
)


# Divide and conquer MDS
groups = sample(x = 5, size = nrow(x), replace = TRUE)

mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)



procrustes_result =  smacof::Procrustes(
  X = mds_classical, 
  Y = mds_divide_conquer$mds
)

df_classical = as.data.frame(mds_classical)
df_classical$type = "classical"
df_classical$label = row.names(x)

df_divide_conquer = as.data.frame(procrustes_result$Yhat)
df_divide_conquer$type = 'divide_conquer'
df_divide_conquer$label = row.names(x)

df_all = rbind(
  df_classical,
  df_divide_conquer
) %>% 
  arrange(
    as.numeric(label)
  )

# Plot
df_all[1:20,] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Euclidean distance
distance_classical_divide = diag(
  rdist(
    mds_classical,
    procrustes_result$Yhat
  )
)

# Plot the error
ggplot(
  data.frame(
    error = distance_classical_divide
  ), 
  aes(error)
) +
  geom_density()

summary(distance_classical_divide)



# 02 Iris data set ----
x = iris[, -5]

dist_matrix = daisy(
  x = x,
  metric = "euclidean"
)


mds_classical = stats::cmdscale(
  d = dist_matrix,
  k = 2
)


# Divide and conquer MDS
groups = sample(x = 5, size = nrow(x), replace = TRUE)

mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)



procrustes_result =  smacof::Procrustes(
  X = mds_classical, 
  Y = mds_divide_conquer$mds
)

df_classical = as.data.frame(mds_classical)
df_classical$type = "classical"
df_classical$label = row.names(x)

df_divide_conquer = as.data.frame(procrustes_result$Yhat)
df_divide_conquer$type = 'divide_conquer'
df_divide_conquer$label = row.names(x)

df_all = rbind(
  df_classical,
  df_divide_conquer
) %>% 
  arrange(
    as.numeric(label)
  )

# Plot
df_all%>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Euclidean distance
distance_classical_divide = diag(
  rdist(
    mds_classical,
    procrustes_result$Yhat
  )
)

# Plot the error
ggplot(
  data.frame(
    error = distance_classical_divide
  ), 
  aes(error)
) +
  geom_density()

summary(distance_classical_divide)





# 03 MNIST data set ----



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
row.names(eastD)
attr(eastD, "Labels") <- abbreviate(attr(eastD, "Labels"))
fit.east <- mds(eastD, type = "ordinal")
westD <- sim2diss(EW_eng$west)
attr(westD, "Labels") <- abbreviate(attr(westD, "Labels"))
fit.west <- mds(westD, type = "ordinal", init = torgerson(eastD))

fit.proc <- Procrustes(fit.east$conf, fit.west$conf)
fit.proc

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
row.names(fit.proc)
row.names()