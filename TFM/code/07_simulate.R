# This one will do 18 replicas
source("tools/load_libraries.R")
source("tools/classical_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/divide_conquer_mds.R")
source("tools/gower_interpolation_mds.R")
source("tools/simulator.R")
source("tools/compute_accuracy.R")

threshold_main_dimensions = 0.9

simulation_id = 40000
initial_scenario_id = 40000
total_replicas = 18


df = expand.grid(
  scenario_id = list(NULL),
  sample_size = list(10^6),
  data_dimension = list(10,100),
  main_dimensions_vector = list(NULL, 15, c(15,15), c(15,10), c(15,15,15,15)),
  l = list(500),
  k = list(3),
  metric = list("euclidean"),
  compute_divide_conquer_mds = list(TRUE),
  compute_fast_mds = list(TRUE),
  compute_gower_mds = list(TRUE),
  compute_classical_mds = list(TRUE),
  max_sample_size_classical = 3000,
  n_eigenvalues = 6,
  n_cols_procrustes_noise = 5, # When there is noise, use 5 columns to do the procrustes
  split_procrustes = list(TRUE), # This is because when the matrix is too big and procrustes is perform to get the correlation matrix, it explodes. With that, procrustes is calcualted in a for
  n_max_procrustes = list(5000)
)

df$scenario_id = initial_scenario_id:(initial_scenario_id + nrow(df)-1)


if(FALSE){
  View(df)
  df = df[c(9,10, 11, 17, 18, 19), ]
  i_row = 2
}


nrows_df = nrow(df)

for(i_replica in 1:total_replicas){
  message("------------------------------------------------------------------------")
  
  message(
    paste0(
      "------- Working on replica ", i_replica, " out of ", total_replicas, "-------"
    )
  )
  
  i_row_ini = 1
  for(i_row in i_row_ini:nrows_df){
    df_filter = df[i_row, ]
    
    # Security control
    if( df_filter$sample_size[[1]] > df_filter$max_sample_size_classical[[1]] ) {
      df_filter$compute_classical_mds[[1]] = FALSE
    }
    
    message(
      paste0(
        "---- Working on scenario ", i_row," out of ", nrows_df, " ----"
      )
    )
    
    n_primary_dimensions = length(df_filter$main_dimensions_vector[[1]])
    n_secondary_dimensions = df_filter$data_dimension[[1]] - n_primary_dimensions
    
    list_results_i = do.magic(
      sample_size = df_filter$sample_size[[1]],
      data_dimension = df_filter$data_dimension[[1]], 
      n_primary_dimensions = n_primary_dimensions,
      
      main_dimensions_vector = df_filter$main_dimensions_vector[[1]],
      l = df_filter$l[[1]],
      k = df_filter$k[[1]],
      metric = df_filter$metric[[1]],
      compute_divide_conquer_mds = df_filter$compute_divide_conquer_mds[[1]],
      compute_fast_mds = df_filter$compute_fast_mds[[1]],
      compute_classical_mds = df_filter$compute_classical_mds[[1]],
      compute_gower_mds = df_filter$compute_gower_mds[[1]],
      max_sample_size_classical = df_filter$max_sample_size_classical[[1]],
      threshold_main_dimensions = threshold_main_dimensions,
      n_eigenvalues = df_filter$n_eigenvalues[[1]],
      n_cols_procrustes_noise =  df_filter$n_cols_procrustes_noise[[1]],
      split_procrustes = df_filter$split_procrustes[[1]],
      n_max_procrustes = df_filter$n_max_procrustes[[1]]
    )
    
    
    
    exists_dominant_dimesion = FALSE
    if( length( df_filter$main_dimensions_vector[[1]] ) > 1){
      exists_dominant_dimesion = length( unique(df_filter$main_dimensions_vector[[1]]) ) > 1
    }
    

    df_summary_i = data.frame(
      scenario_id = df_filter$scenario_id,
      simulation_id = simulation_id,
      
      sample_size = list_results_i$sample_size,
      n_dimensions = list_results_i$data_dimension,
      n_primary_dimensions = n_primary_dimensions,
      value_primary_dimensions = NA,
      n_secondary_dimensions = n_secondary_dimensions,
      
      exists_dominant_dimesion = exists_dominant_dimesion,
      
      # Output for divide and conquer
      n_dimensions_divide_conquer = list_results_i$divide_conquer_n_dimensions,
      elapsed_time_divide_conquer = list_results_i$divide_conquer_elapsed_time,
      eig_subsample_divide_conquer = NA,
      corr_matrix_divide_conquer = NA,
      
      # Output for fast
      n_dimensions_fast = list_results_i$fast_n_dimensions,
      elapsed_time_fast = list_results_i$fast_elapsed_time,
      eig_subsample_fast = NA,
      corr_matrix_fast = NA,
      
      # Output form gower
      gower_n_dimensions = list_results_i$gower_n_dimensions,
      elapsed_time_gower = list_results_i$gower_elapsed_time,
      eig_subsample_gower = NA,
      corr_matrix_gower = NA,
      
      # Output for classical
      n_dimensions_classical = list_results_i$classical_n_dimensions,
      elapsed_time_classical = list_results_i$classical_elapsed_time,
      eig_subsample_classical = NA,
      corr_matrix_classical = NA
    )
    df_summary_i$value_primary_dimensions = list(list_results_i$main_dimensions_vector)
    
    df_summary_i$eig_subsample_divide_conquer = list(list_results_i$divide_conquer_eig_subsample)
    df_summary_i$corr_matrix_divide_conquer = list(list_results_i$divide_conquer_corr_matrix)

    df_summary_i$eig_subsample_fast = list(list_results_i$fast_eig_subsample)
    df_summary_i$corr_matrix_fast = list(list_results_i$fast_corr_matrix)
    
    df_summary_i$eig_subsample_gower = list(list_results_i$gower_eig_subsample)
    df_summary_i$corr_matrix_gower = list(list_results_i$gower_corr_matrix)

    df_summary_i$eig_subsample_classical = list(list_results_i$classical_eig_subsample)
    df_summary_i$corr_matrix_classical = list(list_results_i$classical_corr_matrix)
    
    
    if(i_row  == 1 && i_replica == 1){
      df_summary = df_summary_i
    }else{
      df_summary = rbind(
        df_summary,
        df_summary_i
      )
    }
    
    simulation_id = simulation_id + 1
    save(df_summary, file = "df_summary.RData")  
    save(df, file = "df.RData")
  }
}

if(FALSE){
  View(df_summary)
  
  df_summary %>% 
    group_by(
      scenario_id
    ) %>% 
    summarise(
      n()
    ) %>% 
    View
}