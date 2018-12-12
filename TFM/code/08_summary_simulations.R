source("tools/load_libraries.R")
load("ws_gr_2018_12_10_19_09_59.Rproj")


View(df_summary_all)

df_summary_all %>% 
  arrange(
    sample_size_divide_conquer_fast,
    n_dimensions,
    n_primary_dimensions
  ) %>% 
  View


# Dimesionality
df_summary_all %>% 
  arrange(
    sample_size_divide_conquer_fast,
    n_dimensions,
    n_primary_dimensions
  ) %>% 
  select(
    sample_size_divide_conquer_fast,
    sample_size_classical,
    n_dimensions,
    n_primary_dimensions,
    n_secondary_dimensions,
    n_dimensions_classical,
    n_dimensions_divide,
    n_dimensions_fast
  ) %>% 
  View

# Timing
df_timing_classical = df_summary_all %>% 
  mutate(
    sample_size = sample_size_classical
  ) %>% 
  group_by(
    sample_size
  ) %>% 
  summarise(
    elapsed_time_classical_mean = mean(elapsed_time_classical)
  )
  
df_timing_new_approaches = df_summary_all %>% 
  mutate(
    sample_size = sample_size_divide_conquer_fast
  ) %>% 
  group_by(
    sample_size
  ) %>% 
  summarise(
    elapsed_time_divide_conquer_mean = mean(elapsed_time_divide_conquer),
    elapsed_time_fast_mean = mean(elapsed_time_fast)
  )


