source("tools/divide_change_partitions.R")
library(tidyverse)

n_sim = 500
sample_size = 1000
n_cols = 10
sd_sim = 10
l = 100
s = n_cols + 1
k = n_cols
increment_vector = seq(0, 0.5, 0.01)
df = data.frame(simulation_id=rep(NA, n_sim*length(increment_vector)), 
                num_partitions=rep(NA, n_sim*length(increment_vector)), 
                elapsed_time=rep(NA, n_sim*length(increment_vector)),
                default_num_partitions=rep(NA, n_sim*length(increment_vector)),
                decreasement=rep(NA, n_sim*length(increment_vector)))

i_iteration = 1
for(increment in increment_vector){
  print(paste0("Starting with increment ", increment))
  for(i_sim in 1:n_sim){
    x = matrix(data=rnorm(n=sample_size*n_cols, sd=sd_sim), nrow=sample_size, ncol=n_cols)
    init_time = proc.time()
    results = divide_conquer_mds(x=x, l=l, s=s, k=k, increment=increment)
    elapsed_time = (proc.time() - init_time)[3]
    
    df$simulation_id[i_iteration] = stringi::stri_rand_strings(n=1, length=15)
    df$num_partitions[i_iteration] = results$p
    df$elapsed_time[i_iteration] = elapsed_time
    df$default_num_partitions[i_iteration] = results$default_num_partitions
    df$decreasement[i_iteration] = increment
    
    i_iteration = i_iteration + 1
  }
}

df %>% 
  group_by(decreasement) %>% 
  summarise(time = median(elapsed_time)) %>% 
  ggplot(aes(x=decreasement, y=time)) +
  geom_line()


  group_by(num_partitions) %>% 
  summarise(n = n()) %>% View
