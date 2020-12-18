source("tools/load_libraries.R")
source("tools/mds_methods.R")
source("tools/procrustes.R")

validate_scenarios <- function(df, what){
  if("keys" %in% what){
    expected_keys = c("sample_size", "n_cols", "distribution_parameters", "mu", "sd", "processed_at")
    key_names = colnames(df)
    for(ek in expected_keys){
      if(!ek %in% key_names) {stop(paste0(ek, " should be a key in scenarios"))}
    }
  }
  
  if("content" %in% what){
    total_scenarios = nrow(df)
    for(i in 1:total_scenarios){
      current_scenario = df[i, ]
      if(current_scenario$n_cols < length(current_scenario$mu)) {stop("There is one scenario with more mean values than columns")}
      if(current_scenario$n_cols < length(current_scenario$sd)) {stop("There is one scenario with more sd values than columns")}
    }
  }
}


#'@title Describe scenarios
#'@description Lists all the scenarios to be simulated
#'@param sample_size List of all sample sizes
#'@param distribution_parameters List of all (mean, sd) 
#'@return Returns a data frame that contains all the scenarios
generate_df_scenarios <- function(scenarios){
  
  df = expand.grid(scenarios)
  df$id = 1:nrow(df)
  df = df[, ncol(df):1]
  df$mu = NA
  df$sd = NA
  df$n_main_dimensions = NA
  df$processed_at = as.POSIXct(NA)
  
  validate_scenarios(df=df, what="keys")
  
  total_scenarios = dim(df)[1]
  for(i in 1:total_scenarios){
    current_scenario = df[i, ]
    n_cols = current_scenario$n_cols
    distribution_parameters = current_scenario$distribution_parameters[[1]]
    mu = unlist(distribution_parameters$mu)
    sd = unlist(distribution_parameters$sd)
    n_main_dimensions = 0
    
    if(is.null(mu)){
      mu = rep(0, times=n_cols)
    }else if(length(mu) < n_cols){
      mu = c(mu, rep(0, times=n_cols-length(mu)))
    }
    
    if(is.null(sd)){
      sd = rep(1, times=n_cols)
    }else if(length(sd) < n_cols){
      n_main_dimensions = length(sd)
      sd = c(sd, rep(1, times=n_cols-length(sd)))
    }
    
    df$mu[i] = list(mu)
    df$sd[i] = list(sd)
    df$n_main_dimensions[i] = n_main_dimensions
  }
  validate_scenarios(df=df, what="content")
  return(df)
}


#'@title Generate data
#'@description Generate a multi-dimensional Normal Distribution
#'@param sample_size List of all sample sizes
#'@param distribution_parameters List of all (mean, sd) 
#'@return Returns a matrix with sample_size rows and distribution_parameters columns
generate_data <- function(scenario){

  total_columns = scenario$n_cols
  sample_size = scenario$sample_size
  mu = unlist(scenario$mu)
  sd = unlist(scenario$sd)
  
  return(mapply(rnorm, n=sample_size, mean=mu, sd=sd))
}


create_correlation_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created){
    load(file_path)
  }else{
    df_correlation = data.frame(scenario_id=numeric(0), num_sim=numeric(0), method_name=character(0),
                                n_main_dimensions=numeric(0), correlation_vector=numeric(0))
  }
  
  #In case a simulation breaks in the middle we delete all the simulations
  processed_scenarios = df_scenarios$id[!is.na(df_scenarios$processed_at)]
  ids_prcoessed = which(df_correlation$scenario_id %in% processed_scenarios)
  df_correlation = df_correlation[ids_prcoessed, ]
  save(df_correlation, file=file_path)
  assign("df_correlation", df_correlation, envir=.GlobalEnv)
  
}


create_time_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
  }else{
    df_time = data.frame(scenario_id=numeric(0), num_sim=numeric(0), method_name=character(0), elapsed_time=numeric(0))
  }
  
  #In case a simulation breaks in the middle we delete all the simulations
  processed_scenarios = df_scenarios$id[!is.na(df_scenarios$processed_at)]
  ids_prcoessed = which(df_time$scenario_id %in% processed_scenarios)
  df_time = df_time[ids_prcoessed, ]
  save(df_time, file=file_path)
  assign("df_time", df_time, envir=.GlobalEnv)
  
}


#'@title Create scenarios
#'@description Create a data frame containing all the scenarios
#'@param file_path Path where file should be stored
#'@param sample_size List of all sample sizes
#'@param distribution_parameters List of all (mean, sd) 
#'@return Store and return a data frame with all the scenarios
create_scenarios_file <- function(file_path, scenarios, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
    validate_scenarios(df=df_scenarios, what=c("content", "keys"))
  }else{
    df_scenarios = generate_df_scenarios(scenarios=scenarios)
    save(df_scenarios, file=file_path)
    assign("df_scenarios", df_scenarios, envir=.GlobalEnv)
  }
}


get_correlation_main_dimesions <- function(x, y, num_dimesions){
  
  if(num_dimesions==0){
    return(NA)
  }else{
    
    corr_vector = c()
    x_main = x[, 1:num_dimesions, drop=FALSE]
    y_main = y[, 1:num_dimesions, drop=FALSE]
    
    x_proc = perform_procrustes(x=x_main, target=y_main, matrix_to_transform=x_main, 
                                translation=FALSE, dilation=FALSE)
    
    for(i_dim in 1:num_dimesions){
      current_corr = cor(x_proc[, i_dim], y_main[, i_dim])
      corr_vector = c(corr_vector, current_corr)
    }
    return(corr_vector)
  }
}


update_scenarios_data <- function(file_path, scenarion_id){
  
  df_scenarios$processed_at[df_scenarios$id==scenarion_id] = Sys.time()
  assign("df_scenarios", df_scenarios, envir=.GlobalEnv) 
  save(df_scenarios, file=file_path)
}


update_time_data <- function(file_path, scenario_id, num_sim, method_name, elapsed_time){
  
  temp_df = data.frame(scenario_id=scenario_id, num_sim=num_sim, method_name=method_name, elapsed_time=elapsed_time)
  df_time = rbind(df_time, temp_df)
  assign("df_time", df_time, envir=.GlobalEnv) 
  save(df_time, file=file_path)
}


update_correlation_data <- function(file_path, scenario_id, num_sim, method_name, n_main_dimensions, correlation_vector){
  
  temp_df = data.frame(scenario_id=scenario_id, num_sim=num_sim, method_name=method_name, 
                       n_main_dimensions=n_main_dimensions)
  
  temp_df$correlation_vector = correlation_vector
  
  df_correlation = rbind(df_correlation, temp_df)
  
  assign("df_correlation", df_correlation, envir=.GlobalEnv) 
  save(df_correlation, file=file_path)
  
}


validate_input <- function(list_inputs){
  parameter_names = names(list_inputs)
  for(name in parameter_names){
    values = list_inputs[[name]]
    if(all(is.na(values)) | all(is.null(values))){
      msg = paste0(name, " must not be NULL nor NA")
      stop(msg)
    }
  }
}

get_simulations <-function(
  scenarios, 
  path,
  mds_methods,
  n_simulations,
  overwrite_simulations=FALSE,
  n_sampling_points = NA,
  largest_matrix_efficient_mds = NA,
  num_mds_dimesions = NA,
  verbose = FALSE
){
  
  validate_input(list(scenarios=scenarios, path=path, mds_methods=mds_methods, 
                      n_simulations=n_simulations, largest_matrix_efficient_mds=largest_matrix_efficient_mds))
  
  scenarios_filename = "df_scenarios.RData"
  time_filename = "df_time.RData"
  correlation_filename = "df_correlation.RData"
  
  create_scenarios_file(file_path=file.path(path, scenarios_filename), scenarios=scenarios, overwrite_simulations=overwrite_simulations)
  create_time_file(file_path=file.path(path, time_filename), overwrite_simulations=overwrite_simulations)
  create_correlation_file(file_path=file.path(path, correlation_filename), overwrite_simulations=overwrite_simulations)
  
  methods_names = sapply(as.list(substitute(mds_methods))[-1], deparse)
  total_methods = length(mds_methods)
  df_missing_scenarios = df_scenarios[is.na(df_scenarios$processed_at),]
  total_scenarios = dim(df_missing_scenarios)[1]
  
  if(total_scenarios == 0) {stop("All scenarios are alredy simulated")}
  
  
  for(i_scenario in 1:total_scenarios){
    if(verbose){
      message("--------------------")
      message(paste0("Starting scenario: ", i_scenario))
    }
    i_sim_method = 1
    
    current_scenario = df_scenarios[i_scenario,]
    
    # Set the parameter values for the methods
    s = ifelse(!is.na(n_sampling_points), n_sampling_points, current_scenario$n_main_dimensions + 1)
    k = ifelse(!is.na(num_mds_dimesions), num_mds_dimesions, current_scenario$n_cols)
    l = largest_matrix_efficient_mds
    
    batch_scenario_ids = c()
    batch_num_sims = c()
    batch_method_names = c()
    batch_elapsed_times = c()
    batch_n_main_dimensions = c()
    batch_correlation_vector = list()
    
    for(i_sim in 1:n_simulations){
      if(verbose & (i_sim == 1 | (i_sim)%%10 == 0)){message(paste0("     Starting simulation: ", i_sim))}
      x = generate_data(scenario=current_scenario)
      i_method = 1
      for(method in mds_methods){
        
        starting_time = proc.time()
        result = method(x=x, l=l, s=s, k=k)
        elapsed_time = (proc.time() - starting_time)[3]
        
        batch_scenario_ids = c(batch_scenario_ids, current_scenario$id)
        batch_num_sims = c(batch_num_sims, i_sim)
        batch_method_names = c(batch_method_names, methods_names[i_method])
        batch_elapsed_times = c(batch_elapsed_times, elapsed_time)
        batch_n_main_dimensions = c(batch_n_main_dimensions, current_scenario$n_main_dimensions)
        
        correlation_vector = get_correlation_main_dimesions(x=x, y=result$points, num_dimesions=current_scenario$n_cols)
        batch_correlation_vector[[i_sim_method]] = correlation_vector
        i_method = i_method + 1
        i_sim_method = i_sim_method + 1
        
      }
    }
    
    update_time_data(file_path=file.path(path, time_filename), scenario_id=batch_scenario_ids, 
                     num_sim=batch_num_sims, method_name=batch_method_names, elapsed_time=batch_elapsed_times)
    update_correlation_data(file_path=file.path(path, correlation_filename), scenario_id=batch_scenario_ids,
                            num_sim=batch_num_sims, method_name=batch_method_names, 
                            n_main_dimensions=batch_n_main_dimensions, correlation_vector=batch_correlation_vector)
    update_scenarios_data(file_path=file.path(path, scenarios_filename), scenarion_id=current_scenario$id)
  }
}
