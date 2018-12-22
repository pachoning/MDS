source("tools/load_libraries.R")

correct_ids = function(
  mod,
  col_value
){
  return(col_value%%mod + 30000)
}

is_30000_fixed = FALSE

directory = file.path(
  "data",
  "simulations",
  "1000000"
)

list_files = list.files(directory)
list_summary = grep( pattern = "summary", list_files, value = TRUE)
list_scenarios = setdiff(list_files, list_summary)

for(i_file in list_summary){
  load(file.path(directory, i_file))
  # Para el de 30000 hay que coger a partir de la fila 37
  is_df_30000 = grepl(pattern = "30000", x = i_file)
  if(is_df_30000 == TRUE & is_30000_fixed == FALSE){
    df_summary = df_summary[37:nrow(df_summary), ]
  }
  
  df_summary$scenario_id = correct_ids(
    mod = min(df_summary$scenario_id),
    col_value = df_summary$scenario_id
  )
  
  save(
    df_summary,
    file = file.path(directory, i_file)
  )
    
  rm(df_summary)
}



for(i_file in list_scenarios){
  load(file.path(directory, i_file))
  df$scenario_id = correct_ids(
    mod = min(df$scenario_id),
    col_value = df$scenario_id
  )
  
  save(
    df,
    file = file.path(directory, i_file)
  )
  
  rm(df)
}

