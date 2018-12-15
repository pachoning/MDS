n.dimensions <- function(
  list_eigenvectors,
  threshold_main_dimensions
){
  if( length(list_eigenvectors) <= 1 & 
      (
        is.null(list_eigenvectors) == TRUE || is.na(list_eigenvectors) == TRUE
      )
    ){
    above_threshold = NA
  }else{
    if( is.list(list_eigenvectors) == FALSE ){
      list_eigenvectors = list(list_eigenvectors)
    }

    above_threshold = rapply(
      list_eigenvectors, 
      function(x, th = threshold_main_dimensions) min(which(cumsum(x)/sum(x) > th)), 
      how = "unlist"
    )
    
    above_threshold = max(above_threshold)
  }
  
  return(above_threshold)
}
