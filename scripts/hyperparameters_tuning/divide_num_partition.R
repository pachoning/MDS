source("tools/divide_change_partitions.R")
library(tidyverse)

increment_vector = seq(0, 0.9, 0.01)
n_cols = c(5, 6)
n_sim = 100
sample_size = 1000

sd_sim = 10
l = 100

#s = 2*n_cols
#k = n_cols

df = data.frame(simulation_id=rep(NA, n_sim*length(increment_vector)), 
                num_partitions=rep(NA, n_sim*length(increment_vector)), 
                elapsed_time=rep(NA, n_sim*length(increment_vector)),
                default_num_partitions=rep(NA, n_sim*length(increment_vector)),
                decreasement=rep(NA, n_sim*length(increment_vector)),
                mean_num_points=rep(NA, n_sim*length(increment_vector)))

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
    df$mean_num_points[i_iteration] = results$mean_num_points
    
    i_iteration = i_iteration + 1
  }
}

df = df %>% filter(!is.na(simulation_id))
df$s = s
df$sample_size = sample_size
df$k = k
df$l = l
df$l_lower = (1-df$decreasement)*df$l

df_group = df %>% 
  group_by(decreasement) %>% 
  summarise(time = mean(elapsed_time))


df_group %>% 
  ggplot(aes(x=decreasement, y=time)) +
  geom_line() +
  ggtitle('method=divide_conquer; points=1000; cols=10; l=100; k=10;')


decr = df_group$decreasement[which.min(df_group$time)]

View(df[df$decreasement == decr, ])
head(df[df$decreasement == decr, ])
View(df)
mean(df$elapsed_time[df$decreasement == decr])
decr_4_n_cols = max(df$decreasement[df$mean_num_points >= 4*n_cols & df$mean_num_points < 5*n_cols])
cond = df$decreasement == decr_4_n_cols
View(df[cond, ])
mean(df$elapsed_time[cond])
sum(cond)


df_num_partitions = df
save(df_num_partitions, file="decisions/df_num_partitions.RData")
