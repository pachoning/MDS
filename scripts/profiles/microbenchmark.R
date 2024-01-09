library(microbenchmark)
library(dplyr)

# In the package documentation they do not recommend to use it with large amaunt of code, which is our case:
# This function is only meant for micro-benchmarking small pieces of source code 
# and to compare their relative performance characteristics.

# According to the documentation of the package "o achieved this, the sub-millisecond (supposedly nanosecond)".
# We do not need to go up to the nanosecond. In our experiments, we can see observe differences bertween the methods
# in the seconds level.
times <- 100
df <- data.frame(
  team=rep(c('A', 'B'), each=500),
  points=rnorm(1000, mean=20)
)


# Using microbenchmark ----
results <- microbenchmark(
  aggregate(df$points, list(df$team), FUN=mean),
  df %>% group_by(team) %>% summarise_at(vars(points), list(name = mean)),
  times = times
)

summary_res <- summary(results)
summary_res$mean/1e6

# Using proc.time ----
time_agg <- c()
time_dplyr <- c()
sys_time <- c()

for (i in 1:times) {
  init_time <- proc.time()
  int <- Sys.time()
  res <- aggregate(df$points, list(df$team), FUN=mean)
  end_time <- proc.time()
  et <- Sys.time() - int
  time_agg <- c(time_agg, end_time[3] - init_time[3])
  sys_time <- c(sys_time, et)
  
  init_time <- proc.time()
  res <- df %>% group_by(team) %>% summarise_at(vars(points), list(name = mean))
  end_time <- proc.time()
  time_dplyr <- c(time_dplyr, end_time[3] - init_time[3])
}

mean(time_agg)
mean(time_dplyr)
