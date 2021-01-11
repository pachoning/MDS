source("tools/simulator.R")

results_folder_name = "experiment_04"
experiment_label = "10^5"

sample_size = c(10*5)
n_cols = c(10, 100)

distribution_parameters = list(list(sd=c(NA)), list(sd=15), list(sd=c(15, 15)), 
                               list(sd=c(15, 10)), list(sd=c(15, 15, 15, 15)))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  experiment_label=experiment_label, 
  scenarios=scenarios,
  path=file.path(getwd(), "data", "experiments", results_folder_name),
  mds_methods_names=c("divide_conquer", "fast", "gower"),
  n_simulations=100,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=200,
  largest_matrix_efficient_procrustes=5000,
  num_mds_dimesions=NA,
  verbose=TRUE
)
