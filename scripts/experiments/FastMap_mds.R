library(Rdimtools)
library(bigmds)

n_rows <- 1000000
n_col <- 10

data <- matrix(data = rnorm(n = n_rows * n_col), nrow = n_rows)
init_time <- proc.time()
fastmap_mds <- do.fastmap(data)
end_time <- proc.time()
end_time[3] - init_time[3]

init_time_2 <- proc.time()
divide_conquer_mds <- divide_conquer_mds(data, l = 100, r = 2, c_points = 10)
end_time_2 <- proc.time()
end_time_2[3] - init_time_2[3]


?do.procrustes


zzz <- bigmds::interpolation_mds( data, l = 100, r = 5)
