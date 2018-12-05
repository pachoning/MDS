source("tools/load_libraries.R")
source("tools/classical_mds.R")
source("tools/fast_mds_eigen.R")
source("tools/divide_conquer_mds.R")
source("tools/simulator.R")
source("tools/n_dimensions.R")

threshold_main_dimensions = 0.9
df = expand.grid(
  sample_size = list(500, 1000, 1500),
  data_dimension = list(4,10,100),
  main_dimensions_vector = list(NULL, 5, c(5,5)),
  l = list(1000),
  k = list(3),
  metric = list("euclidean"),
  compute_divide_conquer_mds = list(TRUE),
  compute_fast_mds = list(TRUE),
  compute_classical_mds = list(TRUE),
  sample_size_classical = list(NULL)
)

nrows_df = nrow(df)
for(i_row in 1:nrows_df){
  df_filter = df[i_row, ]
  
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
    sample_size_classical = df_filter$sample_size_classical[[1]]
  )
  
  
  df_summary_i = data.frame(
    sample_size = list_results_i$sample_size,
    
    n_columns = list_results_i$data_dimension,
    
    n_dimensions = ifelse(
      length(list_results_i$main_dimensions_vector) == 0,
      list_results_i$data_dimension,
      length(list_results_i$main_dimensions_vector)
    ),
    
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
    )
  )
  
  if(i_row == 1){
    df_summary = df_summary_i
  }else{
    df_summary = rbind(
      df_summary,
      df_summary_i
    )
  }
  
}

View(df_summary)
