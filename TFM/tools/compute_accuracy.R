compare.methods <- function(
  mds_new_approach,
  mds_classic
){
  
  # Obtaining Procrustes
  procrustes_result =  MCMCpack::procrustes(
    X = mds_classic,
    Xstar = mds_new_approach, #target matrix
    translation = TRUE, 
    dilation = TRUE
  )
  
  
  # Transforming into data frames
  df_new_approach = as.data.frame(mds_new_approach)
  df_classic_mds_transformed = as.data.frame(procrustes_result$X.new)
  
  # Formating the data frame
  total_col = ncol(mds_classic)
  colnames(df_new_approach) = paste0("V", 1:total_col)
  colnames(df_classic_mds_transformed) = paste0("V", 1:total_col)
  
  
  df_classic_mds_transformed$type = "classic"
  df_classic_mds_transformed$label = row.names(x)
  
  
  df_new_approach$type = "new_approach"
  df_new_approach$label = row.names(x)
  
  
  # Concatenating both data frames
  df_all = rbind(
    df_classic_mds_transformed,
    df_new_approach
  ) %>% 
    arrange(
      as.numeric(label)
    )
  
  
  # Euclidean distance
  distance_between_coordinates = diag(
    rdist(
      procrustes_result$X.new,
      mds_new_approach
    )
  )
  
  
  return(
    list(
      mds_classic = mds_classic,
      mds_classic_transformed = procrustes_result$X.new,
      df_both_mds_labels = df_all,
      distance_between_coordinates = distance_between_coordinates
    )
  )
}
