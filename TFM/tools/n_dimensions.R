n.dimensions <- function(
  list_eigenvectors,
  threshold_main_dimensions
){
  if( is.list(list_eigenvectors) == FALSE ){
    list_eigenvectors = list(list_eigenvectors)
  }
  
  list_cumsum = list_eigenvectors %>% map(cumsum)
  total_variance = list_eigenvectors %>% map(sum)
  above_threshold = map2(
    list_cumsum, 
    total_variance,
    ~ .x / .y
    ) %>% 
    map_int(
      function(x, th = threshold_main_dimensions) min(which(x >= th))
    ) %>% 
    mean()
  return(above_threshold)
}
