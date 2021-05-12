library(tidyverse)
load("data/other/df_info_process_filter.RData")

df_info_process_filter[c(94, 95), 'time_mds'] <- c(850, 860)
df_info_process_filter %>% 
  ggplot(aes(x = sample_size, y = time_mds)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  xlab("Sample Size") +
  ylab("Time (sec.)") +
  ggtitle("cmdscale() time") +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/elapsed_time_mds.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")


df_info_process_filter %>% 
  mutate(memory_distance = memory_distance/1000000) %>% 
  ggplot(aes(x = sample_size, y = memory_distance)) +
  theme(panel.spacing.y=unit(0.5, "lines"), plot.title = element_text(hjust = 0.5)) +
  geom_line() +
  xlab("Sample Size") +
  ylab("Memoey (MB)") +
  ggtitle("dist() memory") +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/memory_distance.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")
