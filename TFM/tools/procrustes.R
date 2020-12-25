#'@title Procrustes
#'@description Perfroms Procrustes transformation.
#'@param x The matrix to be transformed.
#'@param target The target matrix.
#'@param k Number of principal coordinates.
#'@param matrix_to_transform The matrix to apply the transformation.
#'@param translation Logical value indicating whether X should be translated.
#'@param dilation logical value indicating whether X should be dilated.
#'@return Returns Procruster to \emph{matrix_to_transform}.
perform_procrustes <- function(x, target, matrix_to_transform, translation, dilation, largest_matrix_efficient_procrustes=5000){
  
  p = 1
  if(nrow(x)>largest_matrix_efficient_procrustes){
    p = ceiling(nrow(x)/largest_matrix_efficient_procrustes) + 2
  }
  
  indexes_group = sort(sample(x=p, size=nrow(x), replace=TRUE))
  
  # Procrustes parameters
  rotation_matrix = matrix(data=0, nrow=ncol(x), ncol=ncol(x))
  translation_matrix = matrix(data=0, nrow=nrow(matrix_to_transform), ncol=ncol(matrix_to_transform))
  cum_dilation = 0
  
  for(i_partition in 1:p){
    
    indexes_current_group = which(indexes_group == i_partition)
    
    procrustes_result =  MCMCpack::procrustes(
      X=x[indexes_current_group, ,drop=FALSE], #The matrix to be transformed
      Xstar=target[indexes_current_group, ,drop=FALSE], # target matrix
      translation=translation, 
      dilation=dilation
    )
  
    rotation_matrix = rotation_matrix + procrustes_result$R
  
    if(translation){
      trans = procrustes_result$tt
    }else{
      trans = matrix(data=0, nrow=ncol(x), ncol=1)
    }
  
    ones_vector = matrix(data=1, nrow=nrow(matrix_to_transform), ncol=1)
    translation_matrix = translation_matrix + ones_vector %*% t(trans)
    
    if(dilation){
      dilation_factor = procrustes_result$s
    }else{
      dilation_factor = 1
    }
    
    cum_dilation = cum_dilation + dilation_factor
  }

  return(cum_dilation/p * matrix_to_transform %*% rotation_matrix/p + translation_matrix/p)
}
