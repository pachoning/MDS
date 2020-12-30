source("tools/mds_methods.R")
source("tools/procrustes.R")

n_rows = 14000
n_cols = 10
x = matrix(rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
var(x[, 9])
dim(x)

s = n_cols + 5
k = n_cols
l = 100
divide_results = divide_conquer_mds(x=x,l=l, s=s, k=k)
divide_proc = perform_procrustes(x=divide_results$points, target=x, matrix_to_transform=divide_results$points, 
                                 translation=FALSE, dilation=FALSE)

cor(divide_proc[,1], x[, 1])

fas_results = fast_mds(x=x,l=l, s=s, k=k)
fast_proc = perform_procrustes(x=fas_results$points, target=x, matrix_to_transform=fas_results$points, 
                                 translation=FALSE, dilation=FALSE)
cor(fast_proc[,1], x[, 1])
