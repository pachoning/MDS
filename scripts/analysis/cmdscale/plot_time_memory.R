library(tidyverse)
load("data/cmdscale/df_info_process_filter.RData")

df_info_process_filter[c(94, 95), 'time_mds'] <- c(850, 860)
df_info_process_filter %>% 
  ggplot(aes(x = sample_size, y = time_mds)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  xlab("Sample Size") +
  ylab("Time (sec.)") +
  ggtitle("cmdscale() time") +
  ggsave(file.path(getwd(), "images", "plot_time_memory.png"),
         dpi = 300, dev = 'png', height = 16, width = 20, units = "cm")

df_info_process_filter %>% 
  mutate(memory_distance = memory_distance/1000000) %>% 
  ggplot(aes(x = sample_size, y = memory_distance)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  xlab("Sample Size") +
  ylab("Memoey (MB)") +
  ggtitle("dist() memory") +
  ggsave(file.path(getwd(), "images", "memory_distance.png"),
         dpi = 300, dev = 'png', height = 16, width = 20, units = "cm")
