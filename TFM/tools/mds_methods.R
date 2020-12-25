#source("tools/load_libraries.R")
#source("tools/procrustes.R")


#'@title Classical MDS
#'@description Perfroms classical MDS with eigenvalues.
#'@param x Data matrix.
#'@param k Number of principal coordinates.
#'@return Returns MDS as well as eigenvalues.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
classical_mds <- function(x, k, return_distance_matrix=FALSE){
  
  mds = list()
  dist_matrix = dist(x=x)
  mds_result = cmdscale(d=dist_matrix, k=k,eig=TRUE)
  
  mds$points = mds_result$points
  mds$eigen = mds_result$eig
  
  if(return_distance_matrix) {mds$distance = as.matrix(dist_matrix)}
  
  return(mds)
}


#'@title Number of partitions
#'@descriptionReturn the number of partitions.
#'@param n Number of rows.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 1 + estimated datsa dimension.
#'@param k Number of principal coordinates.
#'@return Returns p (number of partitions).
get_partitions_for_mds <- function(n, l, s, k){
  
  p = ceiling(l/s)
  num_obs_group = ceiling(n/p)
  
  if(num_obs_group < k+2){stop("Too many columns and too few observations to perform Fast MDS")}
  
  while(p>0 & num_obs_group < k + 2){
    p = p -1
    num_obs_group = ceiling(n/p)
  }
  
  partition = sort(rep(x=1:p, length.out=n, each=num_obs_group))
  last_partition_value = max(partition)
  last_group_idx = which(partition == last_partition_value)
  
  if(length(last_group_idx) <= k){
    partition[last_group_idx] = last_partition_value-1
  }
  
  return(partition)
}


#'@title Fast MDS
#'@description Perfroms MDS based on Tynia, Jinze, Leonard, and Wei (2006).
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 1 + estimated datsa dimension.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on fast MDS algorithm.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
fast_mds <- function(x,l,s,k,largest_matrix_efficient_procrustes=5000){

  has_row_names = !is.null(row.names(x))
  if(!has_row_names){
    row.names(x) = 1:nrow(x)
  }
  
  #If possible to run classical MDS on the whole matrix, run it
  if(nrow(x)<=l){
    mds = classical_mds(x=x, k=k)
    mds$eigen = mds$eigen/length(mds$eigen)
    
    if(!has_row_names){
      row.names(x) = NULL
      row.names(mds) = NULL
    }
    
    return(mds)
    
    # Otherwise, call it recursively
  }else{
    index_partition = get_partitions_for_mds(n=nrow(x), l=l, s=s, k=k)
    p = length(unique(index_partition))
    points = list()
    eigen = c()
    min_len = Inf
    sampling_points = list()
    
    # For each partition, compute fast MDS
    for(i in 1:p){
      indexes_partition = which(index_partition==i)
      x_partition = x[indexes_partition, ,drop=FALSE]
      mds_partition = fast_mds(x=x_partition, l=l, s=s, k=k)
      points[[i]] = mds_partition$points
      row.names(points[[i]]) = row.names(x_partition)
      sampling_points[[i]] = sample(x=row.names(x_partition), size=s, replace=FALSE)
      
      if(i==1){
        min_len = length(mds_partition$eigen)
        eigen = mds_partition$eigen
      }else{
        min_len = pmin(min_len, length(mds_partition$eigen))
        eigen = eigen[1:min_len] + mds_partition$eigen[1:min_len]
      }
      
    }
    
    # Perform the mean for the eigenvalues
    eigen = eigen/p
    
    # Get M_align by getting the sampled points
    ind = unlist(sampling_points)
    x_M = x[ind, ,drop=FALSE]
    row.names(x_M) = row.names(x[ind, ,drop=FALSE])
    mds_M = classical_mds(x=x_M, k=k)
    mds_M = mds_M$points
    row.names(mds_M) = row.names(x_M)
    
    #Align the solutions
    for(i in 1:p){
      mds_i = points[[i]]
      sampling_points_i = sampling_points[[i]] 
      mds_M_sampling = mds_M[sampling_points_i, ,drop=FALSE]
      mds_i_sampling = mds_i[sampling_points_i, ,drop=FALSE]
      
      mds_aligned_i = perform_procrustes(x=mds_i_sampling, target=mds_M_sampling, 
                                         matrix_to_transform=mds_i, 
                                         translation=FALSE, dilation=FALSE,
                                         largest_matrix_efficient_procrustes=largest_matrix_efficient_procrustes)
      if(i==1){
        mds_stitched = mds_aligned_i
      }else{
        mds_stitched = rbind(mds_stitched, mds_aligned_i)
      }
    } 
    
    if(!has_row_names){
      row.names(x) = NULL
      row.names(mds_stitched) = NULL
    }
    
    return(list(points=mds_stitched, eigen=eigen))
  }
}


get_partitions_for_divide_conquer <- function(n, l, s, k){
  p = ceiling(n/(l-(s+2)))
  index_partition = sort(rep(x=1:p, length.out=n, each=ceiling(n/p)))
  
  if(mean(table(index_partition)) < k+2) stop("Too many columns and too few observations to perform Divide and Conquer MDS")
  
  n_iter = 1
  while((sum(index_partition == p) < k +2) & (n_iter <= 10)){
    p = p + 1
    index_partition = sort(rep(x=1:p, length.out=n, each=ceiling(n/p)))
    n_iter = n_iter + 1
  }
  
  return(index_partition)
}


#'@title Divide and Conquer MDS
#'@description Perfroms Divide and Conquer MDS.
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 1 + estimated datsa dimension.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on fast MDS algorithm.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
divide_conquer_mds <- function(x,l,s,k,largest_matrix_efficient_procrustes=5000){
  initial_row_names = row.names(x)
  row.names(x) = 1:nrow(x)
  
  if(nrow(x)<=l){
    mds_to_return = classical_mds(x=x, k=k)
    mds_to_return$eigen = mds_to_return$eigen/length(mds_to_return$eigen)
  }else{
    if(s>l){stop("s cannot be larger than l")}

    index_partition = get_partitions_for_divide_conquer(n=nrow(x), l=l, s=s, k=k)
    p = max(index_partition)
    
    min_len = Inf
    eigen = c()

    # Calculate mds for each partition and take s poits from each subsample
    for (i in 1:p) {
      indexes_current = which(index_partition==i)
      x_current = x[indexes_current, ,drop=FALSE]
      row_names_current = row.names(x_current)
      list_classical_mds = classical_mds(x=x_current, k=k) 
      mds_current = list_classical_mds$points
    
      if(i == 1){
        cum_mds = mds_current
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
    mds_to_return = list(points=cum_mds, eigen=eigen)
  }
  
  row.names(x) = initial_row_names
  row.names(mds_to_return$points) = initial_row_names
  
  return(mds_to_return)
}


#'@title MDS based on Gower Interpolation formula
#'@description Perfroms MDS based on Gower Inpterpolation formula
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on Gower interpolation formula.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
gower_interpolation_mds <- function(x,l,k,...){
  
  nrow_x = nrow(x)
  p = ceiling(nrow_x/l)
  if(p<1) p = 1
  
  if( p>1 ){
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution = sort(sample(x=p, size=nrow_x, replace=TRUE))
    sample_distribution = sort(sample_distribution)
    
    # Get the first group 
    ind_1 = which(sample_distribution==1)
    n_1 = length(ind_1)
    
    # Do MDS with the first group
    submatrix_data = x[ind_1, ,drop=FALSE]
    mds_eig = classical_mds(x=submatrix_data, k=k, return_distance_matrix=TRUE)
    distance_matrix = mds_eig$distance
    
    M = mds_eig$points
    eigen = mds_eig$eig/nrow(M)
    cum_mds = M
    
    # Calculations needed to do Gower interpolation
    delta_matrix = distance_matrix^2 
    In = diag(n_1)
    ones_vector = rep(1, n_1)
    J = In - 1/n_1*ones_vector %*% t(ones_vector)
    G = -1/2 * J %*% delta_matrix %*% t(J) 
    g_vector = diag(G)
    # S = cov(M)
    S = 1/(nrow(M)-1)*t(M) %*% M
    S_inv = solve(S)
    
    # For the rest of the groups, do the interpolation
    for(i_group in 2:p){
      # Filtering the data
      ind_i_group = which(sample_distribution == i_group)
      submatrix_data = x[ind_i_group, ,drop=FALSE]
      
      # A matrix
      distance_matrix_filter = pdist::pdist(
        X = submatrix_data,
        Y = x[ind_1, ,drop=FALSE]
      )
      
      distance_matrix_filter = as.matrix(distance_matrix_filter)
      A = distance_matrix_filter^2
      ones_vector = rep(1, length(ind_i_group))
      MDS_i_group = 1/(2*n_1)*(ones_vector %*%t(g_vector) - A) %*% M %*% S_inv
      cum_mds = rbind(
        cum_mds,
        MDS_i_group
      )
    }
  }else{
    # It is possible to run MDS directly
    mds_eig = classical_mds(x=x, k=k, return_distance_matrix=TRUE)
    distance_matrix = mds_eig$distance
    
    cum_mds = mds_eig$points
    eigen = mds_eig$eig/nrow_x
  }
  
  return(list(points=cum_mds, eigen=eigen))
}
