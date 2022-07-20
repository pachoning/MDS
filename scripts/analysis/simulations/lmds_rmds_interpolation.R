library(tidyverse)
load('data/experiments/experiment_lmds_emmanuel_10/df_time.RData')
load('data/experiments/experiment_lmds_emmanuel_10/df_scenarios.RData')

df_time$Algorithm <- df_time$method_name
df_time$method_name <- NULL
df_time$Algorithm <- factor(df_time$Algorithm,
                              levels = c('lmds', 'emmanuel', 'gower'),
                              labels = c('LMDS', 'RMDS', 'Interpolation')
                              )
df_data <- df_time %>% 
  left_join(
    df_scenarios %>% select(id, sample_size), 
    by = c("scenario_id"="id")
  ) %>% 
  mutate(
    log_elapsed = log10(elapsed_time),
    log_sample_size = log10(sample_size)
  )

df_data %>%
  group_by(sample_size, Algorithm) %>% 
  summarise(time = mean(elapsed_time)) %>% 
  pivot_wider(names_from = Algorithm, values_from=time)

df_data %>%
  group_by(log_sample_size, Algorithm) %>% 
  summarise(log_time = mean(log_elapsed)) %>% 
  pivot_wider(names_from = Algorithm, values_from=log_time)

df_data %>%
  group_by(sample_size, Algorithm) %>% 
  summarise(time = mean(elapsed_time)) %>% 
  ggplot(aes(x = sample_size, y = time, group = Algorithm, color = Algorithm)) +
  geom_line() + 
  geom_point() +
  xlab("Sample size") +
  ylab("Elapsed time (sec.)")


df_data %>%
  group_by(log_sample_size, Algorithm) %>% 
  summarise(log_time = mean(log_elapsed)) %>% 
  ggplot(aes(x = log_sample_size, y = log_time, group = Algorithm, color = Algorithm)) +
  geom_line() + 
  geom_point() + 
  xlab("Log sample size") +
  ylab("Log elapsed time (sec.)")

