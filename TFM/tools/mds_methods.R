source("tools/load_libraries.R")
source("tools/procrustes.R")


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
get_number_partitions <- function(n, l, s, k){
  p = ceiling(l/s)
  reminder = n%%p
  quotient = floor(n/p)
  while((reminder!=0 & reminder<=k) | quotient<=k){
    p = p-1
    reminder = n%%p
    quotient = floor(n/p)
  }
  return(p)
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
fast_mds <- function(x,l,s,k){
  
  has_row_names = !is.null(row.names(x))
  if(!has_row_names){
    row.names(x) = 1:nrow(x)
  }
  
  #If possible to run classical MDS on the whole matrix, run it
  if(nrow(x)<=l){
    
    mds = classical_mds(x=x, k=k)
    
    if(!has_row_names){
      row.names(x) = NULL
      row.names(mds_stitched) = NULL
    }
    
    return(mds)
  
  # Otherwise, call it recursively
  }else{
    p = get_number_partitions(n=nrow(x), l=l, s=s, k=k)
    index_partition = rep(x=1:p, length.out=nrow(x), each=floor(nrow(x)/p))
    points = list()
    eigen = list()
    sampling_points = list()
    
    # For each partition, compute fast MDS
    for(i in 1:p){
      indexes_partition = which(index_partition==i)
      x_partition = x[indexes_partition, ]
      mds_partition = fast_mds(x=x_partition, l=l, s=s, k=k)
      points[[i]] = mds_partition$points
      row.names(points[[i]]) = row.names(x_partition)
      eigen[[i]] = mds_partition$eigen
      sampling_points[[i]] = sample(x=row.names(x_partition), size=s, replace=FALSE)
    }
    
    # Get M_align by getting the sampled points
    ind = unlist(sampling_points)
    x_M = x[ind, ]
    row.names(x_M) = row.names(x[ind, ])
    mds_M = classical_mds(x=x_M, k=k)
    mds_M = mds_M$points
    row.names(mds_M) = row.names(x_M)
    
    #Align the solutions
    for(i in 1:p){
      mds_i = points[[i]]
      sampling_points_i = sampling_points[[i]] 
      mds_M_sampling = mds_M[sampling_points_i,]
      mds_i_sampling = mds_i[sampling_points_i, ]
    
      mds_aligned_i = perform_procrustes(x=mds_i_sampling, target=mds_M_sampling, 
                                         matrix_to_transform=mds_i, 
                                         translation=FALSE, dilation=FALSE)
      if(i==1){
        mds_stitched = mds_aligned_i
        eigen_stitched = eigen[[i]]
      }else{
        mds_stitched = rbind(mds_stitched, mds_aligned_i)
        eigen_stitched = c(eigen_stitched, eigen[[i]])
      }
    } 
    
    if(!has_row_names){
      row.names(x) = NULL
      row.names(mds_stitched) = NULL
    }
    
    return(list(points=mds_stitched, eigen=eigen_stitched))
  }
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
divide_conquer_mds <- function(x,l,s,k){
  
  initial_row_names = row.names(x)
  row.names(x) = 1:nrow(x)
  
  if(nrow(x)<=l){
    mds_to_return = classical_mds(x=x, k=k)
  }else{
    
    if(s>l){stop("s cannot be larger than l")}

    p = ceiling(nrow(x)/l)
    index_partition = rep(x=1:p, length.out=nrow(x), each=ceiling(nrow(x)/p))
    eigen = c()

    # Calculate mds for each partition and take s poits from each subsample
    for (i in 1:p) {
      indexes_current = which(index_partition==i)
      x_current = x[indexes_current, ]
      row_names_current = row.names(x_current)
      list_classical_mds = classical_mds(x=x_current, k=k) 
      mds_current = list_classical_mds$points
      row.names(mds_current) = row.names(x_current)
      eigen = c(eigen, list_classical_mds$eigen)
    
      if(i == 1){
        cum_mds = mds_current
      }else{
        indexes_previous = which(index_partition==(i-1))
        row_names_previous = row.names(x)[indexes_previous]
        rn_subsample_previous = sample(x=row_names_previous, size=s, replace=FALSE)
        
        list_mds_both = classical_mds(x=x[c(rn_subsample_previous, row_names_current), ], k=k)
        mds_both = list_mds_both$points
        row.names(mds_both) = c(rn_subsample_previous, row_names_current)
        mds_both_previous = mds_both[rn_subsample_previous, ]
        mds_both_current = mds_both[row_names_current, ]
        cum_mds_previous = cum_mds[rn_subsample_previous, ]
        mds_current_aligned = perform_procrustes(
          x=mds_both_previous, target=cum_mds_previous, matrix_to_transform=mds_both_current, 
          translation=FALSE, dilation=FALSE
        )
        row.names(mds_current_aligned) = row_names_current
        cum_mds = rbind(cum_mds, mds_current_aligned)
      }
    }
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
gower_interpolation_mds <- function(
  x,
  l,
  k,
  ...
){
  
  nrow_x = nrow(x)
  p = ceiling(nrow_x/l)
  if(p<1) p = 1
  
  if( p>1 ){
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution = sample(x=p, size=nrow_x, replace=TRUE)
    sample_distribution = sort(sample_distribution)
    
    # Get the first group 
    ind_1 = which(sample_distribution==1)
    n_1 = length(ind_1)
    
    # Do MDS with the first group
    submatrix_data = x[ind_1, ]
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
      submatrix_data = x[ind_i_group, ]
      
      
      # A matrix
      distance_matrix_filter = pdist::pdist(
        X = submatrix_data,
        Y = x[ind_1, ]
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
  
  return(
    list(
      points = cum_mds,
      eig = eigen
    )
  )
}
