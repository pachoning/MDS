source("final_methods/fast_mds.R")
source("final_methods/interpolation_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/landmark_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")
source("tools/procrustes.R")
source("tools/new_simulator.R")

initial_time <- lubridate::now("UTC")

simulator(
  methods = list(
    divide_conquer_mds = divide_conquer_mds,
    interpolation_mds = interpolation_mds,
    fast_mds = fast_mds,
    landmark_mds = landmark_mds,
    pivot_mds = pivot_mds,
    reduced_mds = reduced_mds
  ),
  n_samples = list(5000, 6000),
  variances = list(
    c(rep.int(15, 2), rep.int(1, 8)),
    c(rep.int(15, 4), rep.int(1, 6))
  ),
  l_global = list(150, 200),
  n_sim = 4,
  factor_n_points = 5,
  path = file.path("data", "fake")
)

message(paste0("Initial time: ", initial_time, " (UTC)"))
message(paste0("end time: ", lubridate::now("UTC"), " (UTC)"))