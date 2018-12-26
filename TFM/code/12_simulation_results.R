source("tools/load_libraries.R")
library(reshape2)

this_directory = rstudioapi::getActiveDocumentContext()$path
main_directory = dirname(dirname(this_directory))

data_directory = file.path(
  main_directory,
  "data"
)

input_data = file.path(
  data_directory,
  "results",
  "df_simulations.RData"
)

# Load data
load(input_data)

nrow(df_simulations)
df_common = df_simulations[, 1:8]

if(FALSE){
  
  df_simulations %>% 
    group_by(
      scenario_id,
      sample_size,
      n_dimensions,
      n_primary_dimensions,
      exists_dominant_dimesion
    ) %>% 
    summarise(
      n(),
      mean(elapsed_time_divide_conquer),
      mean(elapsed_time_fast),
      mean(elapsed_time_gower)
    ) %>% 
    View
  
}

# Eigenvalues of divide conquer 
df = map_df(
  df_simulations$eig_subsample_divide_conquer, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4],
    val5 = x[5],
    val6 = x[6]
  )
)

df_eigen_divide_conquer = cbind(
  df_common,
  df
)
rm(df)





# Correlation of divide conquer
df = map_df(
  df_simulations$corr_matrix_divide_conquer, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4]
  )
)


df_corr_divide_conquer = cbind(
  df_common,
  df
)
rm(df)



# Eigenvalues of fast
df = map_df(
  df_simulations$eig_subsample_fast, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4],
    val5 = x[5],
    val6 = x[6]
  )
)

df_eigen_fast = cbind(
  df_common,
  df
)
rm(df)



# Correlation fast
df = map_df(
  df_simulations$corr_matrix_fast, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4]
  )
)


df_corr_fast = cbind(
  df_common,
  df
)
rm(df)


# Eigenvalues of Gower
df = map_df(
  df_simulations$eig_subsample_gower, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4],
    val5 = x[5],
    val6 = x[6]
  )
)

df_eigen_gower = cbind(
  df_common,
  df
)
rm(df)

# Correlation Gower
df = map_df(
  df_simulations$corr_matrix_gower, 
  function(x) data.frame(
    val1 = x[1],
    val2 = x[2],
    val3 = x[3],
    val4 = x[4]
  )
)


df_corr_gower = cbind(
  df_common,
  df
)
rm(df)


# Boxplot for correlation divide and conquer
df_melt_corr_divide_conquer = melt(
  df_corr_divide_conquer,
  id = colnames(df_corr_divide_conquer)[1:8]
)


ggplot(
  df_melt_corr_divide_conquer, 
  aes(x=scenario_id, y=value, fill=variable)
) +
  geom_boxplot() +
  facet_wrap(~scenario_id, scale="free")
