#'@title Procrustes
#'@description Perfroms Procrustes transformation.
#'@param x The matrix to be transformed.
#'@param target The target matrix.
#'@param k Number of principal coordinates.
#'@param matrix_to_transform The matrix to apply the transformation.
#'@param translation Logical value indicating whether X should be translated.
#'@param dilation logical value indicating whether X should be dilated.
#'@return Returns Procruster to \emph{matrix_to_transform}.
perform_procrustes <- function(x, target, matrix_to_transform, translation, dilation){
  
  procrustes_result =  MCMCpack::procrustes(
    X=x, #The matrix to be transformed
    Xstar=target, # target matrix
    translation=translation, 
    dilation=dilation
  )
  
  rotation_matrix = procrustes_result$R
  
  if(translation){
    trans = procrustes_result$tt
  }else{
    trans = matrix(data=0, nrow=ncol(x), ncol=1)
  }
  
  ones_vector = matrix(data=1, nrow=nrow(matrix_to_transform), ncol=1)
  translation_matrix = ones_vector %*% t(trans)
  
  if(dilation){
    dilation = procrustes_result$s
  }else{
    dilation = 1
  }

  return(dilation * matrix_to_transform %*% rotation_matrix + translation_matrix)
}
