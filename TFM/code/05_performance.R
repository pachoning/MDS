# 01 Load libraries ----
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_MDS.R")





# 02 Load data set ----


# 03 Params ----

# 04 Performance fast MDS algorithm ----
first_number_groups = 20
last_number_groups = 200
number_coordinates = 2
i_iteration = 1

total_number_groups = last_number_groups - first_number_groups + 1


df_performance = data.frame(
  n_obs = rep(NA, total_number_groups),
  computed_time = rep(NA, total_number_groups),
  max_observations_per_group = rep(NA, total_number_groups),
  mean_observations_per_group = rep(NA, total_number_groups),
  
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
  
  # Observations per group
  df_observations_per_group = df_split %>% 
    group_by(group_member) %>% 
    summarise(
      total_observations = n()
    ) %>% 
    ungroup(
      
    ) %>% 
    summarise(
      min_observations_per_group = min(total_observations),
      max_observations_per_group = max(total_observations),
      mean_observations_per_group = mean(total_observations)
    )
  
  
  # Storing the performance informatiob
  df_performance$number_groups[i_iteration] = i_group
  
  df_performance$min_observations_per_group[i_iteration] = df_observations_per_group$min_observations_per_group
  df_performance$max_observations_per_group[i_iteration] = df_observations_per_group$max_observations_per_group
  df_performance$mean_observations_per_group[i_iteration] = df_observations_per_group$mean_observations_per_group
  
  df_performance$elapsed_time_divide_conquer[i_iteration] = elapsed_time_divide_conquer
  
  df_performance$max_error_procrustes[i_iteration] = max(procrustes_error)
  df_performance$min_error_procrustes[i_iteration] = min(procrustes_error) 
  df_performance$mean_error_procrustes[i_iteration] = mean(procrustes_error)   
  df_performance$median_error_procrustes[i_iteration] = median(procrustes_error)   
    
  
  i_iteration = i_iteration + 1
}

save_directory = file.path(
  getwd(),
  "data",
  "Bike-Sharing-Dataset",
  "df_performance.RData"
)
save(df_performance, file = save_directory)

View(df_performance)
