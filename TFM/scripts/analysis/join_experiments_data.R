library(tidyverse)
library(stringi)

data_folder = file.path(getwd(), 'data')
excluded_experiments = 'experiment_03'

all_files = list.files(data_folder)
all_experiments = all_files[which(stringi::stri_detect_fixed(str=all_files, pattern="experiment_"))]
experiments_to_retrieve = setdiff(all_experiments, excluded_experiments)
total_experiments_to_retrieve = length(experiments_to_retrieve)

df_scenarios_full = c()
df_time_full = c()
df_correlation_full = c()

i_experiment = 1
while(i_experiment <= total_experiments_to_retrieve){
  current_experiment = file.path(data_folder, experiments_to_retrieve[i_experiment])
  experiment_files = list.files(current_experiment)
  for(file in experiment_files){
    if(stringi::stri_detect_regex(str=file, pattern="df_scenarios.RData")){
      load(file.path(current_experiment, file))
      df_scenarios_full = rbind(df_scenarios_full, df_scenarios)
      rm(df_scenarios)
    }else if(stringi::stri_detect_regex(str=file, pattern="df_time.RData")){
      load(file.path(current_experiment, file))
      df_time_full = rbind(df_time_full, df_time)
      rm(df_time)
      
    }else if(stringi::stri_detect_regex(str=file, pattern="df_correlation.RData")){
      load(file.path(current_experiment, file))
      df_correlation_full = rbind(df_correlation_full, df_correlation)
      rm(df_correlation)
    }
  }
  
  i_experiment = i_experiment + 1
  
}
  

# Validations
setdiff(df_scenarios_full$id, unique(df_time_full$scenario_id))
setdiff(unique(df_time_full$scenario_id), df_scenarios_full$id)
setdiff(df_scenarios_full$id, unique(df_correlation_full$scenario_id))
setdiff(unique(df_correlation_full$scenario_id), df_scenarios_full$id)
setdiff(unique(df_time_full$scenario_id), unique(df_correlation_full$scenario_id))
setdiff(unique(df_correlation_full$scenario_id), unique(df_time_full$scenario_id))

save(df_scenarios_full, file=file.path(data_folder, 'df_scenarios_full.RData'))
save(df_time_full, file=file.path(data_folder, 'df_time_full.RData'))
save(df_correlation_full, file=file.path(data_folder, 'df_correlation_full.RData'))
