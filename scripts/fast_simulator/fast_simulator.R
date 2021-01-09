library(tidyverse)
source("tools/mds_methods.R")
source("tools/procrustes.R")

n_rows = 50000
n_cols = 10
n_sim = 50

s = 2*n_cols
k = n_cols
l = 100

df = data.frame(method_name=character(0), 
                elapsed_time=numeric(0))

for(i_sim in 1:n_sim){
  if(i_sim%%10==0 | i_sim==1){message(paste0("Simulation number: ", i_sim))}
  x = matrix(rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
  starting_time = proc.time()
  divide_results = divide_conquer_mds(x=x, l=l, s=s, k=k)
  elapsed_time_divide = (proc.time() - starting_time)[3]
  
  starting_time = proc.time()
  fast_results = fast_mds(x=x,l=l, s=s, k=k)
  elapsed_time_fast = (proc.time() - starting_time)[3]
  
  df_temp = data.frame(method_name=c("divide", "fast"), elapsed_time=c(elapsed_time_divide, elapsed_time_fast))
  df  = rbind(df,df_temp)
  
}

df %>% 
  group_by(method_name) %>% 
  summarise(time=median(elapsed_time))


### 1-0.4
#method_name  time
#<chr>       <dbl>
#1 divide      0.751
#2 fast        0.607