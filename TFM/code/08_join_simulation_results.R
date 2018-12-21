source("tools/load_libraries.R")
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


View(df_simulations)

df_simulations %>% 
  group_by(
    sample_size
  ) %>% 
  summarise(
    n()
  ) 

save(
  df_simulations,
  file = file.path(
    output_directory,
    "df_simulations.RData"
  )
)
