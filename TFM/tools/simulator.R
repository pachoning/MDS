validate_scenarios <- function(scenarios){
  expected_keys = c("sample_size", "n_cols", "distribution_parameters")
  key_names = names(scenarios)
  for(ek in expected_keys){
    if(!ek %in% key_names) {stop(paste0(ek, "should be a key in scenarios"))}
  }
}

#'@title Describe scenarios
#'@description Lists all the scenarios to be simulated
#'@param sample_size List of all sample sizes
#'@param distribution_parameters List of all (mean, sd) 
#'@return Returns a data frame that contains all the scenarios
describe_scenarios <- function(scenarios){

  df = expand.grid(scenarios)
  df$id = 1:nrow(df)
  df = df[, ncol(df):1]
  df$processed_at = NA
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
  distribution_parameters = scenario$distribution_parameters[[1]]
  mu = distribution_parameters$mu
  sd = distribution_parameters$sd
  
  if(length(mu)<total_columns){
    mu_0 = rep(0, times=total_columns-length(mu))
    mu = c(mu, mu_0)
  }
  
  if(length(sd)<total_columns){
    sd_1 = rep(1, times=total_columns-length(sd))
    sd = c(sd, sd_1)
  }
  
  
  return(mapply(rnorm, n=sample_size, mean=mu, sd=sd))
}


#'@title Create scenarios
#'@description Create a data frame containing all the scenarios
#'@param file_path Path where file should be stored
#'@param sample_size List of all sample sizes
#'@param distribution_parameters List of all (mean, sd) 
#'@return Store and return a data frame with all the scenarios
create_scenarios_file <- function(file_path, scenarios){
  
  if(file.exists(file_path)){
    load(file_path)
  }else{
    df_scenarios = describe_scenarios(scenarios=scenarios)
    save(df_scenarios, file=file_path)
    
  }
  
  filter = is.na(df_scenarios$processed_at)
  return(df_scenarios[filter, ])
  
}


get_simulations <-function(
  scenarios, 
  path
){
  
  validate_scenarios(scenarios= scenarios)
  scenarios_filename = "df_scenarios.RData"
  df_scenarios = create_scenarios_file(file_path=file.path(path, scenarios_filename), scenarios=scenarios)
  total_scenarios = dim(df_scenarios)[1]
  if(total_scenarios == 0) {stop("All scenarios are alredy simulated")}
  m = list()
  for(i in 1:total_scenarios){
    
    current_scenario = df_scenarios[i,]
    x = generate_data(
      scenario=current_scenario
    )
    
    m[[i]] = x
  
  }
  return(m)
  
}


sample_size = c(100, 1000)
n_cols = c(10, 100)
distribution_parameters = list(list(mu=5, sd=c(10, 15)), list(mu=7, sd=c(15, 15)), list(sd=10))
scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

sss = get_simulations(
  scenarios=scenarios,
  path='./data'
)
