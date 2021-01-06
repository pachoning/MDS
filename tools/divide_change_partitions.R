source("tools/procrustes.R")
classical_mds <- function(x, k, return_distance_matrix=FALSE){
  
  mds = list()
  dist_matrix = dist(x=x)
  mds_result = cmdscale(d=dist_matrix, k=k,eig=TRUE)
  
  mds$points = mds_result$points
  mds$eigen = mds_result$eig
  
  if(return_distance_matrix) {mds$distance = as.matrix(dist_matrix)}
  
  return(mds)
}


get_partitions_for_divide_conquer <- function(n, l, s, k, increment){
  
  p = ceiling(n/l)
  default_num_partitions = p
  
  l_lower = (1-increment)*l
  min_sample_size = max(k+2, s)
  
  p = ceiling(n/l_lower)

  index_partition = sort(rep(x=1:p, length.out=n, each=ceiling(n/p)))
  p = max(index_partition)

  if(mean(table(index_partition)) < min_sample_size) stop("Too many columns and too few observations to perform Divide and Conquer MDS")
  
  n_iter = 1
  while((min(table(index_partition)) < min_sample_size) & (n_iter <= 10)){
    p = p + 1
    index_partition = sort(rep(x=1:p, length.out=n, each=ceiling(n/p)))
    n_iter = n_iter + 1
  }
  
  if(min(table(index_partition)) < min_sample_size){
    stop("Partitions suffer from lack of data")
  }
  
  return(list(index_partition=index_partition, default_num_partitions=default_num_partitions))
}


divide_conquer_mds <- function(x, l, s, k, increment, largest_matrix_efficient_procrustes=5000){
  initial_row_names = row.names(x)
  row.names(x) = 1:nrow(x)
  
  if(nrow(x)<=l){
    mds_to_return = classical_mds(x=x, k=k)
    mds_to_return$eigen = mds_to_return$eigen/length(mds_to_return$eigen)
    mds_to_return$p = 1
    mds_to_return$default_num_partitions = default_num_partitions
  }else{
    if(s>l){stop("s cannot be larger than l")}
    
    partitions = get_partitions_for_divide_conquer(n=nrow(x), l=l, s=s, k=k, increment=increment)
    index_partition = partitions$index_partition
    p = max(index_partition)
    default_num_partitions = partitions$default_num_partitions
    mean_num_points = mean(table(index_partition))
    
    min_len = Inf
    eigen = c()
    
    # Calculate mds for each partition and take s poits from each subsample
    for (i in 1:p) {
      indexes_current = which(index_partition==i)
      x_current = x[indexes_current, ,drop=FALSE]
      row_names_current = row.names(x_current)
      
      if(i == 1){
        list_classical_mds = classical_mds(x=x_current, k=k) 
        cum_mds = list_classical_mds$points
        eigen = list_classical_mds$eigen/length(list_classical_mds$eigen)
        min_len = length(eigen)
      }else{
        
        list_mds_both = classical_mds(x=x[c(rn_subsample_previous, row_names_current), ,drop=FALSE], k=k)
        mds_both = list_mds_both$points
        mds_both_previous = mds_both[rn_subsample_previous, ,drop=FALSE]
        mds_both_current = mds_both[row_names_current, ,drop=FALSE]
        cum_mds_previous = cum_mds[rn_subsample_previous, ,drop=FALSE]
        mds_current_aligned = perform_procrustes(x=mds_both_previous, target=cum_mds_previous, 
                                                 matrix_to_transform=mds_both_current, 
                                                 translation=FALSE, dilation=FALSE,
                                                 largest_matrix_efficient_procrustes=largest_matrix_efficient_procrustes)
        cum_mds = rbind(cum_mds, mds_current_aligned)
        min_len = pmin(min_len, length(list_mds_both$eigen))
        eigen = eigen[1:min_len] + (list_mds_both$eigen[1:min_len]/length(list_mds_both$eigen))
      }
      
      rn_subsample_previous = sample(x=row_names_current, size=s, replace=FALSE)
      
    }
    
    # Perform the mean for the eigenvalues
    eigen = eigen/p
    
    # Divide by the number of observations
    mds_to_return = list(points=cum_mds, eigen=eigen, p=p, default_num_partitions=default_num_partitions, mean_num_points=mean_num_points)
  }
  
  row.names(x) = initial_row_names
  row.names(mds_to_return$points) = initial_row_names
  
  return(mds_to_return)
}
