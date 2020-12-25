source("tools/simulator.R")

data_folder_name = 'experiment_3'
sample_size = c(1000)
n_cols = c(5, 6)

distribution_parameters = list(list(sd=c(NA)), list(sd=c(5)), list(sd=c(5,10)), list(sd=c(5,10, 15)), 
                               list(sd=c(5,10, 15, 20)), list(sd=c(5,10, 15, 25, 30)))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)


get_simulations(
  scenarios=scenarios,
  path=file.path(getwd(), 'data', data_folder_name),
  mds_methods_vector=c(divide_conquer_mds, fast_mds, gower_interpolation_mds),
  n_simulations=3,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=100,
  largest_matrix_efficient_procrustes=5000,
  num_mds_dimesions=NA,
  verbose=TRUE
)
