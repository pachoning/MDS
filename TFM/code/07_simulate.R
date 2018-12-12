source("tools/load_libraries.R")
source("tools/classical_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/divide_conquer_mds.R")
source("tools/simulator.R")
source("tools/n_dimensions.R")
source("tools/compute_accuracy.R")

threshold_main_dimensions = 0.9

initial_simulation_id = floor(1000000*runif(1))
total_replicas = 75


df = expand.grid(
  scenario_id = list(NULL),
  sample_size = list(1000, 3000, 5000, 10000),
  data_dimension = list(4, 10, 100),
  main_dimensions_vector = list(NULL, 15, c(15,10), c(15, 15)),
  l = list(500),
  k = list(3),
  metric = list("euclidean"),
  compute_divide_conquer_mds = list(TRUE),
  compute_fast_mds = list(TRUE),
  compute_classical_mds = list(TRUE),
  max_sample_size_classical = 3000
)

df$scenario_id = initial_simulation_id:(initial_simulation_id + nrow(df)-1)

if(FALSE){
  View(df)
  df = df[1:4, ]
  i_row = 2
}



# list_results <- list()
nrows_df = nrow(df)
for(i_row in 1:nrows_df){
  df_filter = df[i_row, ]
    
  # Security control
  
  if( df_filter$sample_size[[1]] > df_filter$max_sample_size_classical[[1]] ) {
      df_filter$compute_classical_mds[[1]] = FALSE
    }
    
  message(
    paste0(
      "------- Working on iteration ", i_row," out of ", nrows_df, " -------"
    )
  )
    
  for(i_replica in 1:total_replicas){
    paste0(message("Working on replica ", i_replica, " out of ", total_replicas))
    
    list_results_i = do.magic(
      sample_size = df_filter$sample_size[[1]],
      data_dimension = df_filter$data_dimension[[1]], 
      main_dimensions_vector = df_filter$main_dimensions_vector[[1]],
      l = df_filter$l[[1]],
      k = df_filter$k[[1]],
      metric = df_filter$metric[[1]],
      compute_divide_conquer_mds = df_filter$compute_divide_conquer_mds[[1]],
      compute_fast_mds = df_filter$compute_fast_mds[[1]],
      compute_classical_mds = df_filter$compute_classical_mds[[1]],
      max_sample_size_classical = df_filter$max_sample_size_classical[[1]]
    )
      
    list_results_i$threshold = threshold_main_dimensions
      
    exists_dominant_dimesion = FALSE
    if( length( df_filter$main_dimensions_vector[[1]] ) > 1){
      exists_dominant_dimesion = length( unique(df_filter$main_dimensions_vector[[1]]) ) > 1
    }
      
    
    df_summary_i = data.frame(
      scenario_id = df_filter$scenario_id,
      sample_size_divide_conquer_fast = list_results_i$sample_size,
      

      n_dimensions = list_results_i$data_dimension,
      
      n_primary_dimensions = length(list_results_i$main_dimensions_vector),
      
      n_secondary_dimensions = list_results_i$data_dimension - length(list_results_i$main_dimensions_vector),
    
      
      exists_dominant_dimesion = exists_dominant_dimesion,
      
      n_dimensions_classical = n.dimensions(
        list_eigenvectors = list_results_i$classical_eig,
        threshold_main_dimensions = threshold_main_dimensions
      ),
      
      n_dimensions_divide = n.dimensions(
        list_eigenvectors = list_results_i$divide_conquer_eig,
        threshold_main_dimensions = threshold_main_dimensions
      ),
      
      n_dimensions_fast = n.dimensions(
        list_eigenvectors = list_results_i$fast_eig,
        threshold_main_dimensions = threshold_main_dimensions
      ),
      
      
      classical_canonical_corr =list_results_i$classical_canonical_corr,
      divide_conquer_canonical_corr = list_results_i$divide_conquer_canonical_corr,
      fast_canonical_corr = list_results_i$fast_canonical_corr,
      
      classical_corr_data = list_results_i$classical_correlation_data,
      divide_conquer_corr_data = list_results_i$divide_conquer_correlation_data,
      fast_corr_data = list_results_i$fast_correlation_data,
      divide_conquer_corr_classical_mds = list_results_i$divide_conquer_corr_classical_mds,
      fast_corr_mds_classical = list_results_i$fast_corr_mds_classical,
      
      
      elapsed_time_classical = list_results_i$classical_elapsed_time,
      elapsed_time_divide_conquer = list_results_i$divide_conquer_elapsed_time,
      elapsed_time_fast = list_results_i$fast_elapsed_time
    )
    
    initial_simulation_id = initial_simulation_id + 1
    
    # list_results[[i_row]] = list_results_i
    
    if(i_row  == 1 && i_replica == 1){
      df_summary = df_summary_i
    }else{
      df_summary = rbind(
        df_summary,
        df_summary_i
      )
    }
  }
  
  save(df_summary, file = "df_summary.RData")  
  save(df, file = "df.RData")
  
}


View(df_summary)
View(df)
char_time = gsub(
  pattern = "-|:| ",
  replacement = '_',
  x = Sys.time()
)

save.image(
  file = paste0("ws_",char_time, ".Rproj")
)



if(FALSE){
  # This is an example to show to Pedro
  # Align divide and classical
  x = list_results_i$x
  mds_divide_conquer = list_results_i$divide_conquer_points
  mds_classical = list_results_i$classical_points
  
  results_compare_divide_conquer = compare.methods(
    mds_new_approach = mds_divide_conquer,
    mds_classical = mds_classical
  )
  
  # The plot that we always see
  plot(mds_divide_conquer[,1], results_compare_divide_conquer$mds_classical_transformed[,1])
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  plot(mds_divide_conquer[,2], results_compare_divide_conquer$mds_classical_transformed[,2])
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  
  plot(mds_divide_conquer[,3], results_compare_divide_conquer$mds_classical_transformed[,3])
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  plot(mds_divide_conquer[,4], results_compare_divide_conquer$mds_classical_transformed[,4])
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  plot(mds_divide_conquer[,10], results_compare_divide_conquer$mds_classical_transformed[,10])
  abline(a = 0, b = 1, col = 2, lwd = 2)
  
  
  # Canonical correlation between data and divide and conquer
  cc_results = cc(x, mds_divide_conquer)
  cc_results$cor
  
  # Correlation between data and divide and conquer
  cor(x, mds_divide_conquer)
  
  # Correlation between rotated data and divide and conquer
  results_compare_divide_conquer_x = compare.methods(
    mds_new_approach = mds_divide_conquer,
    mds_classical = x
  )
  
  cor(results_compare_divide_conquer_x$mds_classical_transformed, mds_divide_conquer)
  
  # Correlation between classical MDS and divide and conquer
  cor(mds_classical, mds_divide_conquer)
  
  
  # Correlation between rotated classical MDS and divide and conquer
  cor(results_compare_divide_conquer$mds_classical_transformed, mds_divide_conquer)
  
  # So, I will select the following informtion:
  ## Canonical correlation
  ## Correlation between rotated data and new approaches
  ## Correlation between rotated classical MDS and new approaches
  
}



if(FALSE){
  df_summary %>% View
  df_summary %>% colnames()
  df_summary %>% 
    arrange(
      n_dimensions,
      n_primary_dimensions,
      sample_size_divide_conquer_fast,
      exists_dominant_dimesion
    ) %>% 
    View
  
  View(df)
  View(df_summary)
  
  length(list_results)
  list_results[[1]]$x
}


