source("tools/simulator.R")

results_folder_name = "experiment_09"
experiment_label = "fix_bug"

sample_size = c(1000)
n_cols = c(10, 100)

distribution_parameters = list(list(sd=c(NA)), list(sd=2), list(sd=2:3), list(sd=2:4), 
                               list(sd=2:5), list(sd=2:6), list(sd=2:7),list(sd=2:8), list(sd=2:9), 
                               list(sd=2:10), list(sd=2:11))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  experiment_label=experiment_label, 
  scenarios=scenarios,
  path=file.path(getwd(), "data", "experiments", results_folder_name),
  mds_methods_vector=c(divide_conquer_mds, fast_mds, gower_interpolation_mds),
  n_simulations=100,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=200,
  largest_matrix_efficient_procrustes=5000,
  num_mds_dimesions=NA,
  verbose=TRUE
)
