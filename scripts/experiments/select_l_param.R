source("final_methods/fast_mds.R")
source("final_methods/interpolation_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/landmark_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")
source("tools/procrustes.R")
source("tools/new_simulator.R")

simulator(
  methods = list(
    divide_conquer_mds = divide_conquer_mds,
    interpolation_mds = interpolation_mds,
    fast_mds = fast_mds,
    landmark_mds = landmark_mds,
    pivot_mds = pivot_mds,
    reduced_mds = reduced_mds
    
  ),
  n_samples = list(1000000),
  variances = list(
    c(rep.int(15, 10), rep.int(1, 90))
  ),
  l_global = as.list(c(seq(200, 1000, 100), 1500)),
  n_sim = 10,
  path = file.path("data", "l_param")
)

simulator(
  methods = list(
    divide_conquer_mds = divide_conquer_mds,
    interpolation_mds = interpolation_mds,
    fast_mds = fast_mds,
    landmark_mds = landmark_mds,
    pivot_mds = pivot_mds,
    reduced_mds = reduced_mds
    
  ),
  n_samples = list(1000000),
  variances = list(
    c(rep.int(15, 10), rep.int(1, 90))
  ),
  l_global = as.list(c(50, 100, 150, 200, 400)),
  n_sim = 1,
  path = file.path("data", "l_param_small")
)

simulator(
  methods = list(
    divide_conquer_mds = divide_conquer_mds,
    interpolation_mds = interpolation_mds,
    fast_mds = fast_mds,
    landmark_mds = landmark_mds,
    pivot_mds = pivot_mds,
    reduced_mds = reduced_mds
    
  ),
  n_samples = list(1000000),
  variances = list(
    c(rep.int(15, 10), rep.int(1, 90))
  ),
  l_global = as.list(c(50, 100, 150, seq(200, 1000, 100), 1500)),
  n_sim = 20,
  path = file.path("data", "l_all")
)
