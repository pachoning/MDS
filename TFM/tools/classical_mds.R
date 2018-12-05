classical.mds <- function(
  x,
  s,
  metric
){
  
  # Classical MDS
  distance_matrix = daisy(
    x = x,
    metric = metric
  )
  
  cmd_eig = stats::cmdscale(
    d = distance_matrix, 
    k = s,
    eig = TRUE
  )
  
  
  
  return(
    list(
      points = cmd_eig$points,
      eig = cmd_eig$eig 
    )
  )
  
}