# 01 Load libraries ----
library(tidyverse)
source("tools/divide_conquer_mds.R")
source("tools/split_sample.R")

# 02 Load data set ----
input_directory = file.path(getwd(),"data", "Bike-Sharing-Dataset")
input_file = "df_input.RData"


load(file.path(input_directory, input_file))
dim(df_input)

# 03 Performance algorithm ----
first_number_groups = 50
last_number_groups = 52
number_coordinates = 2
i_iteration = 1

total_number_groups = last_number_groups - first_number_groups + 1


df_performance = data.frame(
  number_groups = rep(NA, total_number_groups),
  elapsed_time_divide_conquer = rep(NA, total_number_groups),
  max_error_procrustes = rep(NA, total_number_groups),
  min_error_procrustes = rep(NA, total_number_groups),
  mean_error_procrustes = rep(NA, total_number_groups),
  median_error_procrustes = rep(NA, total_number_groups)
)


for(i_group in first_number_groups:last_number_groups){
  # Spliting the sample
  df_split = split_sample(
    x = df_input,
    is_random_method = TRUE,
    number_groups = i_group,
    variable = NULL,
    return_data_frame = TRUE
  )
  
  
  # Inputs for divide and conquer algorithm
  groups = df_split$group_member
  
  df_select = df_split %>% 
    select(
      -instant,
      -dteday,
      -group_member
    )
  
  # Setting initial time
  ini_proc = proc.time()
  
  # Applying divide and conquer algorithm
  divide_conquer_mds_result = divide_conquer_mds(
    x = df_select,
    groups = groups,
    number_coordinates = number_coordinates,
    metric = "gower"
  )
  
  procrustes_error = divide_conquer_mds_result$error
  
  # Elapsed time calculation for divide and conquer algorithm
  end_proc = proc.time()
  elapsed_time_divide_conquer = as.numeric((end_proc - ini_proc)[3])
  
  # Storing the performance informatiob
  df_performance$number_groups[i_iteration] = i_group
  df_performance$elapsed_time_divide_conquer[i_iteration] = elapsed_time_divide_conquer
  df_performance$max_error_procrustes[i_iteration] = max(procrustes_error)
  df_performance$min_error_procrustes[i_iteration] = min(procrustes_error) 
  df_performance$mean_error_procrustes[i_iteration] = mean(procrustes_error)   
  df_performance$median_error_procrustes[i_iteration] = median(procrustes_error)   
    
  
  i_iteration = i_iteration + 1
}

