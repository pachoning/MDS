source("tools/load_libraries.R")
library(reshape2)
library(car)

# 01 Load data ----
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

# 02 Eigenvalues of divide conquer ----
df = map_df(
  df_simulations$eig_subsample_divide_conquer, 
  function(x) data.frame(
    phi_1 = x[1],
    phi_2 = x[2],
    phi_3 = x[3],
    phi_4 = x[4],
    phi_5 = x[5],
    phi_6 = x[6]
  )
)

df_eigen_divide_conquer = cbind(
  df_common,
  df
)
rm(df)





# 03 Correlation of divide conquer ----
df = map_df(
  df_simulations$corr_matrix_divide_conquer, 
  function(x) data.frame(
    dimension_1 = x[1],
    dimension_2 = x[2],
    dimension_3 = x[3],
    dimension_4 = x[4]
  )
)


df_corr_divide_conquer = cbind(
  df_common,
  df
)
rm(df)



# 04 Eigenvalues of fast ----
df = map_df(
  df_simulations$eig_subsample_fast, 
  function(x) data.frame(
    phi_1 = x[1],
    phi_2 = x[2],
    phi_3 = x[3],
    phi_4 = x[4],
    phi_5 = x[5],
    phi_6 = x[6]
  )
)

df_eigen_fast = cbind(
  df_common,
  df
)
rm(df)



# 05 Correlation fast ----
df = map_df(
  df_simulations$corr_matrix_fast, 
  function(x) data.frame(
    dimension_1 = x[1],
    dimension_2 = x[2],
    dimension_3 = x[3],
    dimension_4 = x[4]
  )
)


df_corr_fast = cbind(
  df_common,
  df
)
rm(df)


# 06 Eigenvalues of Gower ----
df = map_df(
  df_simulations$eig_subsample_gower, 
  function(x) data.frame(
    phi_1 = x[1],
    phi_2 = x[2],
    phi_3 = x[3],
    phi_4 = x[4],
    phi_5 = x[5],
    phi_6 = x[6]
  )
)

df_eigen_gower = cbind(
  df_common,
  df
)
rm(df)

# 07 Correlation Gower ----
df = map_df(
  df_simulations$corr_matrix_gower, 
  function(x) data.frame(
    dimension_1 = x[1],
    dimension_2 = x[2],
    dimension_3 = x[3],
    dimension_4 = x[4]
  )
)


df_corr_gower = cbind(
  df_common,
  df
)
rm(df)


# 09 Boxplot for correlation divide and conquer ----
df_melt_corr_divide_conquer = melt(
  df_corr_divide_conquer,
  id = colnames(df_corr_divide_conquer)[1:8],
  value.name = "correlation"
)


df_melt_corr_divide_conquer %>% 
  filter(
    n_primary_dimensions == 0 
    # !is.na(correlation) & 
      # sample_size == 10^3
  ) %>% 
ggplot(
  aes(x=scenario_id, y = correlation, fill = variable)
) +
  geom_boxplot() +
  facet_wrap(~scenario_id, scale="free")



# 10 Boxplot for fast ----
df_melt_corr_fast = melt(
  df_corr_fast,
  id = colnames(df_corr_fast)[1:8],
  value.name = "correlation"
)



df_melt_corr_fast %>% 
  filter(
    n_primary_dimensions == 0 &
    # sample_size == 10^6  &
      !is.na(correlation) 
  ) %>% 
  ggplot(
    aes(x=scenario_id, y = correlation, fill = variable)
  ) +
  geom_boxplot() +
  facet_wrap(~scenario_id, scale="free")





# 11 Boxplot for correlation gower ----
df_melt_corr_gower = melt(
  df_corr_gower,
  id = colnames(df_corr_gower)[1:8],
  value.name = "correlation"
)

View(df_melt_corr_gower)
View()


df_melt_corr_gower %>% 
  filter(
    n_primary_dimensions > 0 &
      sample_size == 10^6 &
      !is.na(correlation) 
  ) %>% 
  ggplot(
    aes(x=scenario_id, y = correlation, fill = variable)
  ) +
  geom_boxplot() +
  facet_wrap(~scenario_id, scale="free")




# 12 Metrics related to time ----
df_simulations %>% 
  group_by(
    sample_size,
    n_dimensions
  ) %>% 
  summarise(
    mean_divide_conquer = mean(elapsed_time_divide_conquer),
    mean_fast = mean(elapsed_time_fast),
    mean_gower = mean(elapsed_time_gower),
    
    var_divide_conquer = var(elapsed_time_divide_conquer),
    var_fast = var(elapsed_time_fast),
    var_gower = var(elapsed_time_gower)
  ) 

# Anova with two variables: sample size and columns
colnames(df_simulations)
df_simulations_time = df_simulations %>% 
  select(
    scenario_id,
    sample_size,
    n_dimensions,
    n_primary_dimensions,
    exists_dominant_dimesion,
    elapsed_time_divide_conquer,
    elapsed_time_fast,
    elapsed_time_gower
  )


df_time_melt = melt(
  df_simulations_time,
  id = colnames(df_simulations_time)[1:5]
)


df_time_melt$variable = gsub(
  pattern = "elapsed_time_",
  replacement = '',
  x = df_time_melt$variable
)

df_time_melt = df_time_melt %>% 
  rename(
    algorithm = variable
  )


df_time_melt$log_value = log(df_time_melt$value)
df_time_melt$sample_size = as.factor(df_time_melt$sample_size)
df_time_melt$n_dimensions = as.factor(df_time_melt$n_dimensions)
df_time_melt$algorithm = as.factor(df_time_melt$algorithm)
df_time_melt$exists_dominant_dimesion = as.factor(df_time_melt$exists_dominant_dimesion)
df_time_melt$scenario_id = as.factor(df_time_melt$scenario_id)


anova.test = aov(
  log_value ~ algorithm + sample_size + n_dimensions,
  data = df_time_melt
)

summary(anova.test) 

xtable::xtable(
  summary(anova.test) 
)

# Linear model
linearMod <- lm(
  log_value ~ algorithm + sample_size + n_dimensions, 
  data=df_time_melt
)
xtable::xtable(
  summary(linearMod)
)


# Histogram of the time
df_time_melt_filter = df_time_melt %>% 
  filter(
    sample_size == 1000000
  )

df_time_melt_filter %>% 
  ggplot(
  aes(
    value, 
    colour = algorithm
  )
) +
  geom_density() +
  labs(
    x ="elapsed time in seconds split by scenario_id",
    y = "estimated density"
  ) + 
  facet_wrap(~scenario_id, scale="free")

# Anova one factor using scenario_id
anova.test = aov(
  log_value ~ ,
  data = df_time_melt
)



anova.test2 = aov(
  log_value ~ scenario_id + algorithm,
  data = df_time_melt
)

xtable::xtable(summary(anova.test2))

# 13 Eigenvalues for Divide and Conquer ----
df_eigen_divide_conquer %>% 
  filter(
    n_primary_dimensions == 4 &
      exists_dominant_dimesion == FALSE
  ) %>% 
  group_by(
    scenario_id
  ) %>% 
  summarise(
    total = n(),
    
    # 1
    mean_sqrt_phi_1 = mean(sqrt(phi_1)),
    bias_1 = mean_sqrt_phi_1 - 15,
    MSE_1 = sum(15 - sqrt(phi_1))^2/total,
    
    # 2
    mean_sqrt_phi_2 = mean(sqrt(phi_2)),
    bias_2 = mean_sqrt_phi_2 - 15,
    MSE_2 = sum(15 - sqrt(phi_2))^2/total,
    
    
    # 3
    mean_sqrt_phi_3 = mean(sqrt(phi_3)),
    bias_3 = mean_sqrt_phi_3 - 15,
    MSE_3 = sum(15 - sqrt(phi_3))^2/total,
    
    # 4
    mean_sqrt_phi_4 = mean(sqrt(phi_4)),
    bias_4 = mean_sqrt_phi_4 - 15,
    MSE_4 = sum(15 - sqrt(phi_4))^2/total
    
  ) %>% 
  select(
    -total
  )%>% 
  xtable::xtable()
  


# 14 Eigenvalues for Fast ----
df_eigen_fast %>% 
  filter(
    n_primary_dimensions == 4 &
      exists_dominant_dimesion == FALSE
  ) %>% 
  group_by(
    scenario_id
  ) %>% 
  summarise(
    total = n(),
    
    # 1
    mean_sqrt_phi_1 = mean(sqrt(phi_1)),
    bias_1 = mean_sqrt_phi_1 - 15,
    MSE_1 = sum(15 - sqrt(phi_1))^2/total,
    
    # 2
    mean_sqrt_phi_2 = mean(sqrt(phi_2)),
    bias_2 = mean_sqrt_phi_2 - 15,
    MSE_2 = sum(15 - sqrt(phi_2))^2/total,
    
    
    # 3
    mean_sqrt_phi_3 = mean(sqrt(phi_3)),
    bias_3 = mean_sqrt_phi_3 - 15,
    MSE_3 = sum(15 - sqrt(phi_3))^2/total,
    
    # 4
    mean_sqrt_phi_4 = mean(sqrt(phi_4)),
    bias_4 = mean_sqrt_phi_4 - 15,
    MSE_4 = sum(15 - sqrt(phi_4))^2/total
  ) %>% 
  select(
    -total
  )%>% 
  xtable::xtable()


# 15 MDS based on Gower interpolation formula ----
df_eigen_gower %>% 
  filter(
    n_primary_dimensions == 4 &
      exists_dominant_dimesion == FALSE
  ) %>% 
  group_by(
    scenario_id
  ) %>% 
  summarise(
    total = n(),
    
    # 1
    mean_sqrt_phi_1 = mean(sqrt(phi_1)),
    bias_1 = mean_sqrt_phi_1 - 15,
    MSE_1 = sum(15 - sqrt(phi_1))^2/total,
    
    # 2
    mean_sqrt_phi_2 = mean(sqrt(phi_2)),
    bias_2 = mean_sqrt_phi_2 - 15,
    MSE_2 = sum(15 - sqrt(phi_2))^2/total,
    
    
    # 3
    mean_sqrt_phi_3 = mean(sqrt(phi_3)),
    bias_3 = mean_sqrt_phi_3 - 15,
    MSE_3 = sum(15 - sqrt(phi_3))^2/total,
    
    # 4
    mean_sqrt_phi_4 = mean(sqrt(phi_4)),
    bias_4 = mean_sqrt_phi_4 - 15,
    MSE_4 = sum(15 - sqrt(phi_4))^2/total
  ) %>% 
  select(
    -total
  )%>% 
  xtable::xtable()


