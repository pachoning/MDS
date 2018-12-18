source("tools/load_libraries.R")
load("data/simulations/small_size/df_summary.RData")
load("data/simulations/small_size/df.RData")



View(df_summary[10,])
# Time metrics
df_summary %>% 
  group_by(
    sample_size
  ) %>% 
  summarise(
    elapsed_time_divide_conquer = mean(elapsed_time_divide_conquer),
    elapsed_time_fast = mean(elapsed_time_fast),
    elapsed_time_gower = mean(elapsed_time_gower),
    elapsed_time_classical = mean(elapsed_time_classical)
  )
View(df_summary)

# Dimesionality
df1 = df_summary[9,]
View(df1)
x = matrix(rnorm(1000*10), ncol = 10) %*% diag(c(15, rep(1, 9)))
var(x[,1])
sqrt(df1$eig_subsample_divide_conquer[[1]][1]/df1$eig_subsample_divide_conquer[[1]][6])
sqrt(df1$eig_subsample_classical[[1]][1]/df1$eig_subsample_classical[[1]][6])

