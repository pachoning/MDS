library(tidyverse)

get_scenarios <- function(n_samples, variances, include_l, l_param) {
  if (include_l) {
    scenarios <- expand.grid(n_samples, variances, l_param)
    colnames(scenarios) <- c("n_sample", "variance", "l")
  } else{
    scenarios <- expand.grid(n_samples, variances)
    colnames(scenarios) <- c("n_sample", "variance")
  }
  total_scenarios <- nrow(scenarios)
  scenarios$id <- stringi::stri_rand_strings(n = total_scenarios, length = 20)
  scenarios <- scenarios %>% arrange(n_sample)
  return(scenarios)
}

get_sample_size_from_scenario <- function(scenario) {
  return(unlist(scenario$n_sample))
}

get_variance_from_scenario <- function(scenario) {
  return(unlist(scenario$variance))
}

get_n_dominant_dimensions <- function(scenario) {
  variance <- get_variance_from_scenario(scenario)
  return(sum(variance != 1))
}

build_data <- function(scenario) {
  n_sample <- get_sample_size_from_scenario(scenario)
  variance <- get_variance_from_scenario(scenario)
  n_col <- length(variance)
  variance_matrix <- diag(variance)
  sd_matrix <- sqrt(variance_matrix)
  data <- matrix(
    data = rnorm(n = n_sample * n_col),
    nrow = n_sample,
    ncol = n_col
  )
  data <- data %*% sd_matrix
  return(data)
}

get_correlation_coefficients <- function(X, Y) {
  n_col = ncol(X)
  correlation_vector <- rep(NA, n_col)
  for (i in 1:n_col) {
    correlation_vector[i] <- cor(X[, i], Y[, i])
  }
  return(correlation_vector)
}

get_l_value <- function(include_l, scenario, l_methods, method_name) {
  if (include_l) {
    l_value <- unlist(scenario$l)
  } else{
    l_value <- l_methods[[method_name]]
  }
  return(l_value)
}

validate_l_parameter <- function(l_global, l_methods, methods_names) {
  if(!is.null(l_global) & !is.null(l_methods)) {
    stop("'l_global' and 'l_methods' can not be used at the same time")
  }
  
  if (!is.null(l_methods)) {
    l_methods_name <- names(l_methods)
    for (name in methods_names) {
      if (! name %in% l_methods_name) {
        stop(paste0(name, " must be included in 'l_methods'"))
      }
    }
  }
}

create_dir <- function(path) {
  if (!file.exists(path)) {
    dir.create(path)
  } else {
    stop(paste0(path, " alredy exists"))
  }
}

save_scenario_information <- function(scenarios, l_global, l_methods, path) {
  full_path <- file.path(getwd(), path)
  scenario_table <- "scenarios.RData"
  scenarios$l_global <- list(l_global)
  scenarios$l_methods <- list(l_methods)
  file_path <- file.path(full_path, scenario_table)
  save(scenarios, file = file_path)
}

save_times_information <- function(times, path) {
  full_path <- file.path(getwd(), path)
  time_table <- "times.RData"
  file_path <- file.path(full_path, time_table)
  save(times, file = file_path)
}

save_eigenvalues_information <- function(eigenvalues, path) {
  full_path <- file.path(getwd(), path)
  eigenvalue_table <- "eigenvalues.RData"
  file_path <- file.path(full_path, eigenvalue_table)
  save(eigenvalues, file = file_path)
}

save_correlation_information <- function(correlations, path) {
  full_path <- file.path(getwd(), path)
  correlation_table <- "correlations.RData"
  file_path <- file.path(full_path, correlation_table)
  save(correlations, file = file_path)
}

save_l_param_information <- function(l_param, path) {
  full_path <- file.path(getwd(), path)
  l_param_table <- "l_param.RData"
  file_path <- file.path(full_path, l_param_table)
  save(l_param, file = file_path)
}

save_params <- function(params, path) {
  full_path <- file.path(getwd(), path)
  params_table <- "params.RData"
  file_path <- file.path(full_path, params_table)
  save(params, file = file_path)
}

create_time_object <- function(scenario, i_sim, method_name, elapsed_time) {
  times <- data.frame(
    scenario_id = scenario$id,
    simulation_id = i_sim,
    algorithm = method_name,
    elapsed_time = as.numeric(elapsed_time)
  )
  return(times)
}

create_eigen_object <- function(scenario, i_sim, method_name, eigen_vals) {
  eigenvalues <- data.frame(
    scenario_id = scenario$id,
    simulation_id = i_sim,
    algorithm = method_name
  )
  eigenvalues$eigenvalues <- list(eigen_vals)
  return(eigenvalues)
  
}

create_corr_object <- function(scenario, i_sim, method_name, correlation_coeffs) {
  correlations <- data.frame(
    scenario_id = scenario$id,
    simulation_id = i_sim,
    algorithm = method_name
  )
  correlations$correlation_coeffs <- list(correlation_coeffs)
  return(correlations)
}

create_l_param_object <- function(scenario, i_sim, method_name, l_value) {
  l_param_info <- data.frame(
    scenario_id = scenario$id,
    simulation_id = i_sim,
    algorithm = method_name,
    l_value = l_value
  )
  return(l_param_info)
}

create_params_object <- function(scenario, s, tie, factor_n_points) {
  params <- data.frame(
    scenario_id = scenario$id,
    s = s,
    tie = tie,
    factor_n_points = factor_n_points
  )
  return(params)
}

simulator <- function(
    methods,
    n_samples,
    variances,
    n_sim,
    path,
    factor_n_points = NULL,
    l_global = NULL,
    l_methods = NULL,
    s = NULL,
    tie = NULL
){
  
  # Store some initial values
  initial_s <- s
  initial_tie <- tie
  
  # Get the methods
  methods_names <- names(methods)
  total_methods <- length(methods_names)
  
  # Check l parameter
  validate_l_parameter(l_global, l_methods, methods_names)
  
  # Check if l must be taken into account in the simlations
  include_l <- !is.null(l_global)
  
  # Build a matrix containing the scenarios
  scenarios <- get_scenarios(n_samples, variances, include_l, l_global)
  total_scenarios <- nrow(scenarios)
  
  # Store scenarios
  create_dir(path)
  save_scenario_information(scenarios, l_global, l_methods, path)
  
  if (is.null(factor_n_points)) {
    stop("Choose a value for factor_n_points")
  }
  
  # For each scenario
  for (i_scenario in 1:total_scenarios){
    msg_s <- FALSE
    msg_tie <- FALSE
    message(paste0("####### Working on scenario ", i_scenario, " out of ", total_scenarios, " ######"))
    # Take the scenario
    scenario <- scenarios[i_scenario, ]
    
    # Get the number of dominant dimensions
    n_dominant_dimensions <- get_n_dominant_dimensions(scenario)
    
    # For each simulation
    for (i_sim in 1:n_sim) {
      if (i_sim == 1 | i_sim %% 1 == 0 | i_sim == n_sim ) {
        message(
          paste0(
            "\t Working on simulation ", i_sim, " out of ", n_sim, " at ", lubridate::now("UTC"), " (UTC)"
          )
        )
      }
      # Built the data
      data <- build_data(scenario)
      
      # For each method
      for (i_method in 1:total_methods){
        # Take the method
        method_name <- methods_names[i_method]
        # Get l param
        l_value <- get_l_value(include_l, scenario, l_methods, method_name)
        
        # Initialise the parameters of al the methods
        if (is.null(initial_s)) {
          s <- factor_n_points * n_dominant_dimensions
        }
        
        if (l_value/s <= 1) {
          s <- 2 * n_dominant_dimensions
          if (!msg_s) {
            message("\t s is high. Setting to ", s)
            msg_s <- TRUE
          }
        }
        
        if (is.null(initial_tie)) {
          tie <- factor_n_points * n_dominant_dimensions
        }
        if ((l_value - tie) < tie) {
          tie <- 2 * n_dominant_dimensions
          if (!msg_tie) {
            message("\t tie is high. Setting to ", tie) 
            msg_tie <- TRUE
          }
        } 
        
        # Get MDS configuration
        current_method <- methods[[method_name]]
        message(paste0("\t\t Working on method: ", method_name, " using l: ", l_value))
        
        init_time <- proc.time()
        mds <- R.utils::doCall(
          current_method,
          x = data,
          X = data,
          l = l_value,
          s = s,
          tie = tie,
          num_landmarks = l_value,
          pivots = l_value,
          k = n_dominant_dimensions,
          ndim = n_dominant_dimensions,
          n_cores = 1
        )
        end_time <- proc.time()
        points <- mds$points
        
        # Align MDS and the data
        data_selected_columns <- data[, 1:n_dominant_dimensions]
        mds_after_procrustes <- perform_procrustes(
          x = points,
          target = data_selected_columns,
          matrix_to_transform = points
        )
        
        # Create times, eigenvalues and correlation coefficients object
        elapsed_time <- end_time[3] - init_time[3]
        eigen_val <- mds$eigen
        corr_coeffs <- get_correlation_coefficients(X = mds_after_procrustes, Y = data_selected_columns)
        
        time_object <- create_time_object(scenario, i_sim, method_name, elapsed_time)
        eigen_object <- create_eigen_object(scenario, i_sim, method_name, eigen_val)
        corr_oject <- create_corr_object(scenario, i_sim, method_name, corr_coeffs)
        l_param_object <- create_l_param_object(scenario, i_sim, method_name, l_value)
        
        if (i_scenario == 1 & i_sim == 1 & i_method == 1) {
          df_time <- time_object
          df_eigen <- eigen_object
          df_corr <- corr_oject
          df_l_param <- l_param_object
        } else {
          df_time <- rbind(df_time, time_object)
          df_eigen <- rbind(df_eigen, eigen_object)
          df_corr <- rbind(df_corr, corr_oject)
          df_l_param <- rbind(df_l_param, l_param_object)
        }
        
        if (i_sim %% 1 == 0 | i_sim == n_sim) {
          #Store the time
          save_times_information(df_time, path)
          
          # Store the eigenvalues
          save_eigenvalues_information(df_eigen, path)
          
          # Store the correlation coefficients
          save_correlation_information(df_corr, path)
          
          #Store the value of l
          save_l_param_information(df_l_param, path)
        }
        
      }
    }
    # Create params object
    params_object <- create_params_object(scenario, s, tie, factor_n_points)
    if (i_scenario == 1) {
      df_params <- params_object
    } else{
      df_params <- rbind(df_params, params_object)
    }
    # Store the parameters
    save_params(df_params, path)
  }
}
