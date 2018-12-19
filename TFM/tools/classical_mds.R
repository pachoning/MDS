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
  
  eigenvalues_classical =  cmd_eig$eig/nrow(x)
  
  
  
  return(
    list(
      points = cmd_eig$points,
      eig = eigenvalues_classical
    )
  )
  
}