# 01 Load libraries ----
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_MDS.R")


# 02 Load data set ----

#US population
load(file.path("data", "us-population-by-zip-code", "df_pop.RData"))
load(file.path("data", "results", "df_metrics_fast.RData"))

# 03 EDA ----
# Missings data inputation
df_pop_filter = df_pop %>% 
  filter(
    !is.na(minimum_age) & !is.na(maximum_age) & !is.na(gender) & !is.na(zipcode) 
    & !is.na(year) 
  )
df_pop_filter$gender = as.factor(df_pop_filter$gender)
df_pop_filter$zipcode = as.factor(df_pop_filter$zipcode)
dim(df_pop_filter)
head(df_pop_filter)

# 03 Params ----
x_original = df_pop_filter %>% select(-geo_id, -zipcode)
head(x_original)

# For fast MDS
data_set = 'usa_pap'
l = 1000
s = 2
k = 3
metric = 'gower'


# 04 Performance fast MDS algorithm ----

sample_selection = seq(2.5*10^5, 10^6, 5*10^4) 


df_iterative_fast = data.frame(
  data_set = character(length(sample_selection)),
  nrow = numeric(length(sample_selection)),
  ncol = numeric(length(sample_selection)),
  l = numeric(length(sample_selection)),
  s = numeric(length(sample_selection)),
  k = numeric(length(sample_selection)),
  time_computing = numeric(length(sample_selection)),
  stringsAsFactors = FALSE
)


for(i_sam in 1:length(sample_selection)){
  message("Performing iteration: ", i_sam," --------------")
  current_sample_size = sample_selection[i_sam]
  x = x_original %>% slice(1:current_sample_size) 
  
  initial_time = proc.time()
  
  # Fast MDS
  mds_fast = fast_mds(
    x = x,
    n = nrow(x),
    l = l,
    s = s,
    k = k,
    metric = metric
  )
  
  end_time = proc.time()
  elapsed_time = end_time - initial_time
  
  
  df_iterative_fast$data_set[i_sam] = data_set
  df_iterative_fast$nrow[i_sam] = nrow(x)
  df_iterative_fast$ncol[i_sam] = ncol(x)
  df_iterative_fast$l[i_sam] = l
  df_iterative_fast$s[i_sam] = s
  df_iterative_fast$k[i_sam] = k
  df_iterative_fast$time_computing[i_sam] = elapsed_time[3]
  
  df_metrics_fast = rbind(
    df_iterative_fast,
    df_metrics_fast
  )
  
  if( i_sam%%100000 == 0 ){
    save(
      df_metrics_fast,
      file = file.path("data", "results", "df_metrics_fast.RData")
    )
  }
}

