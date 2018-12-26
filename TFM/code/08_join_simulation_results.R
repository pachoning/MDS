source("tools/load_libraries.R")
library(sampling)
# Join all the information
this_directory = rstudioapi::getActiveDocumentContext()$path
main_directory = dirname(dirname(this_directory))

data_directory = file.path(
  main_directory,
  "data"
)

simulations_directory = file.path(
  data_directory,
  "simulations"
)
  
output_directory = file.path(
  data_directory,
  "results"
)

input_list_directories = list.files(simulations_directory)
for(i_directory in input_list_directories){
  message(paste0("Working on directory: ", i_directory))
  current_directory = file.path(
    simulations_directory,
    i_directory
  )
  
  data_lo_be_loaded = grep(
    pattern = "summary",
    list.files(current_directory),
    value = TRUE
  )
  
  for(file_to_be_loaded in data_lo_be_loaded){
    load(
      file.path(
        current_directory,
        file_to_be_loaded
      )
      
    )
    
    if(exists("df_simulations") == FALSE){
      df_simulations = df_summary
    }else{
      df_simulations = rbind(
        df_simulations,
        df_summary
      )
    }
    
    rm(df_summary)
  }
}



  
set.seed(12345)
sss = sampling::strata(
  data = df_simulations,
  stratanames = "scenario_id",
  size = rep(100, length(unique(df_simulations$scenario_id))),
  method = 'srswor'
)
  
df_simulations = df_simulations[sss$ID_unit,]
  
df_simulations = df_simulations %>% 
  arrange(
    scenario_id
  )

View(df_simulations)

nrow(df_simulations)
df_simulations %>% 
  group_by(
    scenario_id,
    sample_size,
    n_dimensions,
    n_primary_dimensions,
    exists_dominant_dimesion
  ) %>% 
  summarise(
    total = n()
  ) %>% 
  View

save(
  df_simulations,
  file = file.path(
    output_directory,
    "df_simulations.RData"
  )
)


if(FALSE){
  df_simulations %>% 
    filter(
      sample_size == 10^6 &
        scenario_id == 30000 
    ) %>% 
    dplyr::select(
      simulation_id,
      scenario_id,
      sample_size,
      n_dimensions,
      n_primary_dimensions,
      value_primary_dimensions,
      n_secondary_dimensions,
      elapsed_time_fast,
      elapsed_time_gower,
      corr_matrix_fast,
      corr_matrix_gower,
      eig_subsample_fast,
      eig_subsample_gower
    ) %>% 
    View
}

