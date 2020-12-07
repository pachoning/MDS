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
classical_mds <-function(x, k){
  
  mds = list()
  dist_matrix = dist(x=x)
  mds_result = cmdscale(d=dist_matrix, k=k,eig=TRUE)
  
  mds$points = mds_result$points
  mds$eigen = mds_result$eig
  
  return(mds)
}


#'@title Number of partitions
#'@descriptionReturn the number of partitions.
#'@param n Number of rows.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 1 + estimated datsa dimension.
#'@param k Number of principal coordinates.
#'@return Returns p (number of partitions).
get_number_partitions <-function(n, l, s, k){
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
fast_mds <-function(x,l,s,k){
  
  has_row_names = !is.null(row.names(x))
  if(!has_row_names){
    row.names(x) = 1:nrow(x)
  }
  
  #If possible to run classical MDS on the whole matrix, run it
  if(nrow(x)<=l){
    mds = classical_mds(x=x, k=k)
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
                                         translation=FALSE, dilation=TRUE)
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

ssss = fast_mds(x=x,l=100,s=10,k=data_dimension)
ssss_proc = perform_procrustes(x=ssss$points, target=x, matrix_to_transform=ssss$points, 
                   translation=FALSE, dilation=FALSE)


cor(ssss_proc[,1], x[,1])

if(TRUE){
  sample_size = 5000
  data_dimension = 5
  main_dimensions_vector = c(2, 2)
  
  x = matrix(
    rnorm(
      n = sample_size*data_dimension
    ),
    ncol = data_dimension,
    nrow = sample_size
  )
  
  dim(x)
  real_data_dimension = length(main_dimensions_vector)
  lambda_vector = rep(1, data_dimension)
  lambda_vector[1:real_data_dimension] = main_dimensions_vector
  
  x = x %*% diag(lambda_vector)
  dim(x)
}
