library(tidyverse)
load("data/cmdscale/df_info_process_filter.RData")

df_info_process_filter[c(94, 95), 'time_mds'] <- c(850, 860)
df_info_process_filter %>% 
  ggplot(aes(x = sample_size, y = time_mds)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  ylab("Time (sec.)") +
  ggsave(file.path(getwd(), "images", "time_cmdscale.png"))

df_info_process_filter %>% 
  mutate(memory_distance = memory_distance/1000000) %>% 
  ggplot(aes(x = sample_size, y = memory_distance)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  ylab("Memoey (MB)") +
  ggsave(file.path(getwd(), "images", "memory_distance.png"))
