source("tools/load_libraries.R")
source("tools/fast_MDS_eigen.R")

determine_depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,determine_depth,thisdepth=thisdepth+1))))    
  }
}



sample_size = 1000
data_dimension = 5
main_dimesions_vector = c(5)
run_classical_mds = TRUE
fast_l = 1000
fast_k = 2
metric = "euclidean"
threshold_main_dimensions = 0.9
compute_fast_MDS = TRUE
compute_divide_conquer_MDS = TRUE
compute_classical_MDS = TRUE
n_obs_classical = 1000





# Build lambda vecotr (dimesionality)
lambda_vector = rep(1, data_dimension)
real_data_dimension = length(main_dimesions_vector)
if( real_data_dimension > 0){
  if(real_data_dimension > data_dimension){
    stop("main_dimesions_vector is greater than data_dimension")
  }
  lambda_vector[1:real_data_dimension] = main_dimesions_vector
} 

# Buil data
x = matrix(
  rnorm(
    n = sample_size*data_dimension
  ),
  ncol = data_dimension,
  nrow = sample_size
)

row.names(x) = as.character(1:sample_size)

x = x %*% diag(lambda_vector)

# Run the classical in case it is needed

# FAST MDS




mds_fast_sol = fast_eigen_mds(
  x = x,
  n = nrow(x),
  l = fast_l,
  s = data_dimension,
  k = fast_k,
  metric = metric
)

mds_fast = mds_fast_sol$points
mds_fast_eigen = mds_fast_sol$eig

depth_list = determine_depth(mds_fast_eigen)
mds_fast_eigen_flatten = mds_fast_eigen

if(depth_list > 1){
  for(i in 1:(depth_list-1)){
    mds_fast_eigen_flatten = unlist(mds_fast_eigen_flatten, recursive = FALSE)
  }
}

length(mds_fast_eigen_flatten)
length(mds_fast_eigen)

list_cumsum = mds_fast_eigen_flatten %>% map(cumsum)
total_variance = mds_fast_eigen_flatten %>% map(sum)
above_threshold = map2(
  list_cumsum, 
  total_variance,
  ~ .x / .y
) %>% 
  map_int(
    function(x, th = threshold_main_dimensions) min(which(x >= th))
  ) 

mean(above_threshold)
barplot(
  prop.table(table(above_threshold)),
  main = paste0("Sample size: ", length(above_threshold))
  )

