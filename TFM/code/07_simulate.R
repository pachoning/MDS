source("tools/load_libraries.R")
source("tools/classical_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/divide_conquer_mds.R")
source("tools/simulator.R")
source("tools/n_dimensions.R")

threshold_main_dimensions = 0.9


df = expand.grid(
  sample_size = list(10^5),
  data_dimension = list(4, 10, 100),
  main_dimensions_vector = list(NULL, 15, c(15,10), c(15, 15)),
  l = list(500),
  k = list(3),
  metric = list("euclidean"),
  compute_divide_conquer_mds = list(TRUE),
  compute_fast_mds = list(TRUE),
  compute_classical_mds = list(TRUE),
  sample_size_classical = list(NULL)
)


if(FALSE){
  View(df)
  df = df[9:10, ]
}


nrows_df = nrow(df)
list_results <- ls()
for(i_row in 1:nrows_df){
  df_filter = df[i_row, ]
  
  # Security control
  if(
    df_filter$sample_size[[1]] >= 5000 & 
    is.null(df_filter$sample_size_classical[[1]]) == TRUE
  ) {
      df_filter$sample_size_classical[[1]] = 3000
  }
  
    message(
      paste0(
        "------- Working on iteration ", i_row," out of ", nrows_df, " -------"
        )
      )
    
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
    
    list_results_i$threshold = threshold_main_dimensions
    
    exists_dominant_dimesion = FALSE
    if( length( df_filter$main_dimensions_vector[[1]] ) > 1){
      exists_dominant_dimesion = length( unique(df_filter$main_dimensions_vector[[1]]) ) > 1
    }
    
  
  df_summary_i = data.frame(
    sample_size_divide_conquer_fast = list_results_i$sample_size,
    
    sample_size_classical = ifelse(
      is.null(list_results_i$sample_size_classical) == TRUE,
      list_results_i$sample_size,
      list_results_i$sample_size_classical
    ),
    
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
    
    elapsed_time_classical = list_results_i$classical_elapsed_time,
    
    elapsed_time_divide_conquer = list_results_i$divide_conquer_elapsed_time,
    
    elapsed_time_fast = list_results_i$fast_elapsed_time
  )
  
  if(i_row == 1){
    df_summary = df_summary_i
    list_results = list_results_i
  }else{
    df_summary = rbind(
      df_summary,
      df_summary_i
    )
    
    list_results = list(
      list_results,
      list_results_i
    )
  }
  
}



char_time = gsub(
  pattern = "-|:| ",
  replacement = '_',
  x = Sys.time()
)

save.image(
  file = paste0("ws_",char_time, ".Rproj")
)


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
}

