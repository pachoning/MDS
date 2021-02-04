source("tools/simulator.R")

results_folder_name = "experiment_99"
experiment_label = "faster_divide_conquer"

sample_size = c(100000)
n_cols = c(100)

distribution_parameters = list(list(sd=c(15, 15, 15, 15)))

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
  largest_matrix_efficient_procrustes = 5000,
  num_mds_dimesions = NA,
  verbose = TRUE
)
