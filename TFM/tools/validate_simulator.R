source("tools/simulator.R")
library(tidyverse)

sample_size = c(5000)
n_cols = c(5)

distribution_parameters = list(list(sd=c(5,10, 15)))
scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  scenarios=scenarios,
  path='./data',
  mds_methods = c(divide_conquer_mds, fast_mds, gower_interpolation_mds),
  n_simulations = 50,
  overwrite_simulations = TRUE,
  n_sampling_points = NA,
  largest_matrix_efficient_mds = 100,
  num_mds_dimesions = NA,
  verbose = TRUE
)

df_time
df_correlation
df_scenarios

df_time %>% group_by(method_name) %>% summarise(mean_time = mean(elapsed_time))
View(df_time)
