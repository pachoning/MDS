source("tools/simulator.R")

results_folder_name = "experiment_89"
experiment_label = "pedro7"

sample_size = 10^5
n_cols = 10

distribution_parameters = list(list(var=rep(15, 4)))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  experiment_label=experiment_label, 
  scenarios=scenarios,
  path=file.path(getwd(), "data", "experiments", results_folder_name),
  mds_methods_names=c("divide_conquer", "fast", "gower"),
  n_simulations=10,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=1000,
  num_mds_dimesions=NA,
  verbose=TRUE
)
