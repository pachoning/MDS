library(tidyverse)

load('data/experiments_processed/df_conversion.RData')
load('data/experiments_processed/df_time_full.RData')
load('data/experiments_processed/df_scenarios_full.RData')


df_scenarios_full$id[duplicated(df_scenarios_full$id)]


df_time_full$Algorithm <- df_time_full$method_name
df_time_full$method_name <- NULL
df_time_full$Algorithm <- factor(df_time_full$Algorithm,
                            levels = c('lmds', 'emmanuel', 'gower'),
                            labels = c('LMDS', 'RMDS', 'Interpolation'))

df_data <- df_time_full %>%
  left_join(
    df_scenarios_full %>% select(id, sample_size), 
    by = c("scenario_id"="id")
  ) %>% 
  mutate(
    log_elapsed = log10(elapsed_time),
    log_sample_size = log10(sample_size),
  )

df_data %>%
  group_by(sample_size, Algorithm) %>%
  summarise(time = mean(elapsed_time), total = n()) %>% 
  pivot_wider(names_from = Algorithm, values_from=time)

df_data %>%
  group_by(sample_size, Algorithm) %>%
  summarise(time = median(elapsed_time), total = n()) %>%
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


#pdf('images/projection_based_comparison.pdf', width = 8, height = 4)
df_data %>%
  group_by(log_sample_size, Algorithm) %>%
  summarise(log_time = mean(log_elapsed)) %>%
  ggplot(aes(x = log_sample_size, y = log_time, group = Algorithm, color = Algorithm)) +
  geom_line() + 
  geom_point() + 
  xlab("Log sample size") +
  ylab("Log elapsed time (sec.)")
#dev.off()

boxplot(log_elapsed ~ Algorithm*sample_size,data = df_data)

df_data %>%
  ggplot(aes(y=elapsed_time, color=Algorithm)) +
  geom_boxplot() +
  facet_wrap(~sample_size, nrow=4, scales = "free")
