source("tools/procrustes.R")

n_rows = 14000
n_cols = 2
x = matrix(rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)


rot_mat = matrix(data=c(cos(45), -sin(45), sin(45), cos(45)), ncol = n_cols, nrow = n_cols)
y = x %*% rot_mat

proc = perform_procrustes(x = x, target = y, matrix_to_transform = x, translation = FALSE, dilation = FALSE)

head(proc)
head(y)
cor(proc[,1], y[,1])
max(abs(proc - y))
