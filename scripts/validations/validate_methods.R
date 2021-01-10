source("tools/mds_methods.R")
source("tools/procrustes.R")

n_rows = 1000
n_cols = 1
x = matrix(rnorm(n_rows*n_cols, sd = 2), nrow = n_rows)
var(x[, 1])
dim(x)

s = 2*n_cols
k = n_cols
l = 100
divide_results = divide_conquer_mds(x=x,l=l, s=s, k=k)
divide_proc = perform_procrustes(x=divide_results$points, target=x, matrix_to_transform=divide_results$points, 
                                 translation=FALSE, dilation=FALSE)

cor(divide_proc[,1], x[, 1])
var(divide_results$points[, 1])

fast_results = fast_mds(x=x,l=l, s=s, k=k)
fast_proc = perform_procrustes(x=fast_results$points, target=x, matrix_to_transform=fast_results$points, 
                                 translation=FALSE, dilation=FALSE)
cor(fast_proc[,1], x[, 1])
