source("tools/mds_methods.R")
source("tools/procrustes.R")
library(tidyverse)

n_sim = 100
sample_size = 5000
n_cols = 15
sd_sim = 10
l = 100
k = n_cols
sampling_points = n_cols:(3*n_cols)

df = data.frame(simulation_id=rep(NA, n_sim*length(sampling_points)), 
                s=rep(NA, n_sim*length(sampling_points)), 
                correlation=rep(NA, n_sim*length(sampling_points)))


get_correlation_main_dimesions <- function(x, y, num_dimesions, largest_matrix_efficient_procrustes){
  
  if(num_dimesions==0){
    return(NA)
  }else if(num_dimesions == 1){
    return(cor(x[, 1], y[, 1]))
  }else{
    
    corr_vector = c()
    x_main = x[, 1:num_dimesions, drop=FALSE]
    y_main = y[, 1:num_dimesions, drop=FALSE]
    
    x_proc = perform_procrustes(x=x_main, target=y_main, matrix_to_transform=x_main, 
                                translation=TRUE, dilation=FALSE,
                                largest_matrix_efficient_procrustes=largest_matrix_efficient_procrustes)
    
    for(i_dim in 1:num_dimesions){
      current_corr = cor(x_proc[, i_dim], y_main[, i_dim])
      corr_vector = c(corr_vector, current_corr)
    }
    return(corr_vector)
  }
}

i_iteration = 1
for(s in sampling_points){
  print(paste0("Starting with s: ", s))
  for(i_sim in 1:n_sim){
    x = matrix(data=rnorm(n=sample_size*n_cols, sd=sd_sim), nrow=sample_size, ncol=n_cols)
    results = divide_conquer_mds(x=x, l=l, s=s, k=k)
    correlation_vector = get_correlation_main_dimesions(x=x, y=results$points, num_dimesions=ncol(x), 
                                                        largest_matrix_efficient_procrustes=5000)
    df$simulation_id[i_iteration] = stringi::stri_rand_strings(n=1, length=15)
    df$s[i_iteration] = s
    df$correlation[i_iteration] = mean(correlation_vector)
    i_iteration = i_iteration + 1
  }
}


df$method = "divide_conquer"
df$sample_size = sample_size
df$n_cols = 15
df$sd_sim = sd_sim
df$l = l
df$k = k

df_group = df %>% 
  group_by(s) %>% 
  summarise(correlation = mean(correlation))


df_group %>% 
  ggplot(aes(x=s, y=correlation)) +
  geom_line() +
  ggtitle('method=divide_conquer; points=5000; cols=15; l=100; k=15')

df_sampling_points = df
save(df_sampling_points, file=file.path('decisions/df_sampling_points.RData'))
