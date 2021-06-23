source("tools/simulator.R")

results_folder_name = "experiment_99"
experiment_label = "using_do_call"

sample_size = c(10^5)
n_cols = c(100)

distribution_parameters = list(list(var=c(15, 15, 15, 15)))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  experiment_label = experiment_label, 
  scenarios = scenarios,
  path = file.path(getwd(), "data", "experiments", results_folder_name),
  mds_methods_names = c("divide_conquer", "fast", "gower"),
  n_simulations = 10,
  overwrite_simulations = FALSE,
  n_sampling_points = NA,
  largest_matrix_efficient_mds = 200,
  num_mds_dimesions = NA,
  verbose = TRUE
)
