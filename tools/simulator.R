source("tools/procrustes.R")

validate_scenarios <- function(df, what){
  if("keys" %in% what){
    expected_keys = c("sample_size", "n_cols", "distribution_parameters", "mu", "var", "processed_at")
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
      if(current_scenario$n_cols < length(current_scenario$var)) {stop("There is one scenario with more var values than columns")}
    }
  }
}

generate_df_scenarios <- function(scenarios, experiment_label){
  
  df = expand.grid(scenarios)
  df$id = stringi::stri_rand_strings(n=nrow(df), length=15)
  df = df[, ncol(df):1]
  df$mu = NA
  df$var = NA
  df$n_main_dimensions = NA
  df$processed_at = as.POSIXct(NA)
  df$computer_id = Sys.info()["nodename"]
  df$experiment_label = experiment_label
  
  validate_scenarios(df=df, what="keys")
  
  total_scenarios = dim(df)[1]
  for(i in 1:total_scenarios){
    
    n_main_dimensions = 0
    current_scenario = df[i, ]
    n_cols = current_scenario$n_cols
    distribution_parameters = current_scenario$distribution_parameters[[1]]

    mu = unlist(distribution_parameters$mu)
    var = unlist(distribution_parameters$var)
    
    if(any(is.na(mu) | is.null(mu))){
      mu = rep(0, times=n_cols)
    }else if(length(mu) < n_cols){
      mu = c(mu, rep(0, times=n_cols-length(mu)))
    }

    
    if(any(is.na(var) | is.null(var))){
      var = rep(1, times=n_cols)
    }else if(length(var) <= n_cols){
      n_main_dimensions = length(var)
      var = c(var, rep(1, times=n_cols-length(var)))
    }
    
    df$mu[i] = list(mu)
    df$var[i] = list(var)
    df$n_main_dimensions[i] = n_main_dimensions
  }
  validate_scenarios(df=df, what="content")
  return(df)
}


generate_data <- function(scenario){

  total_columns = scenario$n_cols
  sample_size = scenario$sample_size
  mu = unlist(scenario$mu)
  var = unlist(scenario$var)
  
  return(mapply(rnorm, mean = mu, n = sample_size) %*% diag(sqrt(var)))
}


create_correlation_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
    warning('Using simulations that are already saved. Current will be ignored')
  }else{
    df_correlation = data.frame(scenario_id=character(0),
                                num_sim=numeric(0),
                                method_name=character(0),
                                n_main_dimensions=numeric(0), 
                                correlation_vector=numeric(0))
  }
  
  #In case a simulation breaks in the middle we delete all the simulations
  save(df_correlation, file=file_path)
  assign("df_correlation", df_correlation, envir=.GlobalEnv)
  
}

create_eigenvalue_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
    warning('Using simulations that are already saved. Current will be ignored')
  }else{
    df_eigenvalue = data.frame(scenario_id=character(0),
                               num_sim=numeric(0),
                               method_name=character(0),
                               eigenvalue_vector=numeric(0))
  }
  
  #In case a simulation breaks in the middle we delete all the simulations
  save(df_eigenvalue, file=file_path)
  assign("df_eigenvalue", df_eigenvalue, envir=.GlobalEnv)
}


create_mds_parameters_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
    warning('Using simulations that are already saved. Current will be ignored')
  }else{
    df_mds_parameters = data.frame(scenario_id = character(0),
                                   s = numeric(0),
                                   k = character(0),
                                   l = numeric(0))
  }
  
  #In case a simulation breaks in the middle we delete all the simulations
  processed_scenarios = df_scenarios$id[!is.na(df_scenarios$processed_at)]
  ids_prcoessed = which(df_mds_parameters$scenario_id %in% processed_scenarios)
  df_mds_parameters = df_mds_parameters[ids_prcoessed, ]
  save(df_mds_parameters, file=file_path)
  assign("df_mds_parameters", df_mds_parameters, envir=.GlobalEnv)
}

create_time_file <- function(file_path, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    load(file_path)
    warning('Using simulations that are already saved. Current will be ignored')
  }else{
    df_time = data.frame(scenario_id = character(0),
                         num_sim = numeric(0),
                         method_name = character(0),
                         elapsed_time = numeric(0))
  }

  save(df_time, file=file_path)
  assign("df_time", df_time, envir=.GlobalEnv)
  
}


create_scenarios_file <- function(file_path, scenarios, experiment_label, overwrite_simulations){
  
  is_file_created = file.exists(file_path)
  if(is_file_created & !overwrite_simulations){
    warning('Using simulations that are already saved. Current will be ignored')
    load(file_path)
    validate_scenarios(df=df_scenarios, what=c("content", "keys"))
  }else{
    folder_path = dirname(dirname(file_path))
    df_scenarios = generate_df_scenarios(scenarios=scenarios, experiment_label=experiment_label)
    save(df_scenarios, file=file_path)
  }
  assign("df_scenarios", df_scenarios, envir=.GlobalEnv)
}


get_correlation_main_dimesions <- function(x, y, num_dimesions){
  
  if(num_dimesions==0){
    return(NA)
  }else if(num_dimesions == 1){
    return(abs(cor(x[, 1], y[, 1])))
  }else{
    
    corr_vector = c()
    x_main = x[, 1:num_dimesions, drop=FALSE]
    y_main = y[, 1:num_dimesions, drop=FALSE]
    
    x_proc = perform_procrustes(x=x_main, target=y_main, matrix_to_transform=x_main, 
                                translation=FALSE)
    
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


update_mds_parameters_data <- function(file_path, scenario_id, s, k, l){
  
  temp_df = data.frame(scenario_id=scenario_id, s=s, k=k, l=l)
  df_mds_parameters = rbind(df_mds_parameters, temp_df)
  
  assign("df_mds_parameters", df_mds_parameters, envir=.GlobalEnv) 
  save(df_mds_parameters, file=file_path)
}


update_time_data <- function(file_path, scenario_id, num_sim, method_name, elapsed_time){
  
  temp_df = data.frame(scenario_id=scenario_id, num_sim=num_sim, method_name=method_name, elapsed_time=elapsed_time)
  df_time = rbind(df_time, temp_df)
  assign("df_time", df_time, envir=.GlobalEnv) 
  save(df_time, file=file_path)
}


update_correlation_data <- function(file_path, scenario_id, num_sim, method_name, n_main_dimensions, correlation_vector){
  
  temp_df = data.frame(scenario_id = scenario_id,
                       num_sim = num_sim,
                       method_name = method_name, 
                       n_main_dimensions = n_main_dimensions)
  
  temp_df$correlation_vector = correlation_vector
  
  df_correlation = rbind(df_correlation, temp_df)
  
  assign("df_correlation", df_correlation, envir=.GlobalEnv) 
  save(df_correlation, file=file_path)
  
}


update_eigenvalue_data <- function(file_path, scenario_id, num_sim, method_name, eigenvalue_vector){
  
  temp_df = data.frame(scenario_id = scenario_id,
                       num_sim = num_sim,
                       method_name = method_name)
  
  temp_df$eigenvalue_vector = eigenvalue_vector
  
  df_eigenvalue = rbind(df_eigenvalue, temp_df)
  
  assign("df_eigenvalue", df_eigenvalue, envir=.GlobalEnv) 
  save(df_eigenvalue, file=file_path)
  
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


get_initial_number_simulations <- function(df_time, scenario_id, n_simulations) {
  initial_simulation = 1
  if(scenario_id %in% df_time$scenario_id) {
    initial_simulation = min(max(df_time$num_sim[df_time$scenario_id == scenario_id]) + 1, n_simulations)
  }
  return(initial_simulation)
}


get_simulations <-function(
  experiment_label,
  scenarios, 
  path,
  algorithms,
  n_simulations,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=NA,
  num_mds_dimesions=NA,
  verbose=FALSE
){
  
  if(!dir.exists(path)) {
    dir.create(path)
  }
  
  validate_input(list(scenarios = scenarios, 
                      path = path, n_simulations = n_simulations,
                      largest_matrix_efficient_mds = largest_matrix_efficient_mds))
  
  input_parameters = as.list(match.call())
  save(input_parameters, file=file.path(path, 'input_parameters.RData'))
  
  scenarios_filename = "df_scenarios.RData"
  time_filename = "df_time.RData"
  correlation_filename = "df_correlation.RData"
  eigenvalue_filename = "df_eigenvalue.RData"
  mds_parameters_filename = "df_mds_parameters.RData"
  
  create_scenarios_file(file_path = file.path(path, scenarios_filename), 
                        scenarios = scenarios, 
                        experiment_label = experiment_label,
                        overwrite_simulations = overwrite_simulations)
  
  create_time_file(file_path = file.path(path, time_filename),
                   overwrite_simulations = overwrite_simulations)
  
  create_correlation_file(file_path = file.path(path, correlation_filename), 
                          overwrite_simulations = overwrite_simulations)
  
  create_eigenvalue_file(file_path = file.path(path, eigenvalue_filename), 
                         overwrite_simulations = overwrite_simulations)
  
  create_mds_parameters_file(file_path = file.path(path, mds_parameters_filename), 
                             overwrite_simulations = overwrite_simulations)
  
  total_methods = length(algorithms)
  df_missing_scenarios = df_scenarios[is.na(df_scenarios$processed_at),]
  total_scenarios = dim(df_missing_scenarios)[1]
  
  if(total_scenarios == 0) {stop("All scenarios are alredy simulated")}
  
  for(i_scenario in 1:total_scenarios){
    if(verbose){
      message("--------------------")
      message(paste0("Starting scenario: ", i_scenario))
    }
    
    current_scenario = df_missing_scenarios[i_scenario, ,drop=FALSE]
    initial_simulation = get_initial_number_simulations(df_time = df_time, scenario_id = current_scenario$id, 
                                                        n_simulations = n_simulations)

    # Set the parameter values for the methods
    s = ifelse(!is.na(n_sampling_points), 
               n_sampling_points, 
               ifelse(current_scenario$n_main_dimensions>0, 
                      2*current_scenario$n_main_dimensions, 
                      min(0.1*current_scenario$sample_size, 10)))
   
    k = ifelse(!is.na(num_mds_dimesions), 
               num_mds_dimesions, 
               pmax(current_scenario$n_main_dimensions, 1))
    
    l = largest_matrix_efficient_mds
    
    for(i_sim in initial_simulation:n_simulations){
      batch_scenario_ids = c()
      batch_num_sims = c()
      batch_method_names = c()
      batch_elapsed_times = c()
      batch_n_main_dimensions = c()
      batch_correlation_vector = list()
      batch_eigenvalue_vector = list()
      
      x = generate_data(scenario=current_scenario)
      n_row_x = nrow(x)
      i_method = 1
      
      if(verbose){
        message(paste0("     Starting simulation: ", i_sim))
      }
      
      i_algoritm = 1
      while (i_algoritm <= total_methods) {
        
        name = names(algorithms)[i_algoritm]
        starting_time = proc.time()
        result = algorithms[[i_algoritm]](x = x, l = l, tie = s, k = k, dist_fn = stats::dist, s = s)
        elapsed_time = (proc.time() - starting_time)[3]
        
        batch_scenario_ids = c(batch_scenario_ids, current_scenario$id)
        batch_num_sims = c(batch_num_sims, i_sim)
        batch_method_names = c(batch_method_names, name)
        batch_elapsed_times = c(batch_elapsed_times, elapsed_time)
        batch_n_main_dimensions = c(batch_n_main_dimensions, current_scenario$n_main_dimensions)

        correlation_vector = get_correlation_main_dimesions(x=x, y=result$points, 
                                                            num_dimesions=k)

        eigenvalue_vector = result$eigen
        batch_correlation_vector[[i_algoritm]] = correlation_vector
        batch_eigenvalue_vector[[i_algoritm]] = eigenvalue_vector
        i_algoritm = i_algoritm + 1
        
      }
      
      update_time_data(file_path = file.path(path, time_filename), 
                       scenario_id = batch_scenario_ids, 
                       num_sim = batch_num_sims, 
                       method_name = batch_method_names, 
                       elapsed_time=batch_elapsed_times)
      
      update_correlation_data(file_path = file.path(path, correlation_filename), 
                              scenario_id = batch_scenario_ids,
                              num_sim = batch_num_sims, 
                              method_name = batch_method_names, 
                              n_main_dimensions = batch_n_main_dimensions, 
                              correlation_vector = batch_correlation_vector)
      
      update_eigenvalue_data(file_path = file.path(path, eigenvalue_filename), 
                             scenario_id = batch_scenario_ids,
                             num_sim = batch_num_sims, 
                             method_name = batch_method_names,
                             eigenvalue_vector = batch_eigenvalue_vector)

    }
    
    update_mds_parameters_data(file_path = file.path(path, mds_parameters_filename), 
                               scenario_id = current_scenario$id, 
                               s = s, 
                               k = k, 
                               l = l)
    
    update_scenarios_data(file_path = file.path(path, scenarios_filename), 
                          scenarion_id = current_scenario$id)
  }
}
