source("tools/load_libraries.R")

compute_coordinates<- function(
  v_postions
){
  x = v_postions%%8
  y = (v_postions - x)/8
  
  y = ifelse(
    x == 0,
    y -1,
    y
  )
  
  x = ifelse(
    x == 0,
    8,
    x
  )
  x = x -1
  
  corrdinates_matrix = matrix(
    data = c(x,y),
    ncol = 2
  )
  
  return(corrdinates_matrix)
  
  
}


calculate_distance <- function(
  vect1,
  vect2,
  number_pieces
){
  max_dist = 7*sqrt(2)
  
  if(
      ( length(vect1) == 0 & length(vect2) != 0) || 
      ( length(vect1) != 0 & length(vect2) == 0)
    ){
    dists <- matrix(max_dist, nrow = number_pieces, ncol = number_pieces)
  } else if( length(vect1) == 0 & length(vect2) == 0 ){
    dists <- matrix(0, nrow = number_pieces, ncol = number_pieces)
  }else{
    mat1 <- data.frame(vect1)
    mat2 <- data.frame(vect2) 
  
    dists <- pdist(mat1, mat2)
    dists = as.matrix(dists)/max_dist
  
  
    if( nrow(dists) < number_pieces ){
      n_rows_to_add = number_pieces - nrow(dists)
      new_matrix_to_append = matrix(1, nrow = n_rows_to_add , ncol = ncol(dists))
      dists = rbind(
        dists,
        new_matrix_to_append
      )
    
    } else if( ncol(dists) < number_pieces ){
      n_cols_to_add = number_pieces - ncol(dists)
      new_matrix_to_append = matrix(1, nrow = nrow(dists) , ncol = n_cols_to_add)
      dists = cbind(
        dists,
        new_matrix_to_append
      )
    }
  }
  
  if( sum( apply(dists, MARGIN = 1, min) ) < sum( apply(dists, MARGIN = 2, min) ) ){
    return( sum( apply(dists, MARGIN = 1, min ) ) )
  } else{
    return( sum( apply(dists, MARGIN = 2, min ) ) )
  }
}





chess_distance <- function(
  df
){
  
  df_pieces <- data.frame(
    piece = c("K", "Q", "R", "B", "N", "P"),
    initial_number = c(1, 1, 2, 2, 2, 8),
    weight = c(11, 9, 5, 3, 3, 1)
  )

  distance_matrix = matrix(
    data = NA,
    nrow = nrow(df),
    ncol = nrow(df)
  )
  
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
      # Take the chess board for game i and j
      dist_ij = 0
      cheass_board_ij = df[c(i, j), ]
      
      # For each piece, compute the distance
      for(i_piece in 1:6){
      
        df_pieces_filter = df_pieces[i_piece, ]
        # For each colour, compute the distance
        for(i_colour in c("White", "Black")){
        
          ind_i = which( cheass_board_ij[i,] == paste0(df_pieces_filter$piece,i_colour)  ) 
          ind_j = which( cheass_board_ij[j,] == paste0(df_pieces_filter$piece,i_colour)  ) 
        
          coordinates_i = compute_coordinates(ind_i)
          coordinates_j = compute_coordinates(ind_j)
        
          distance_pieces = calculate_distance(
            vect1 = ind_i, 
            vect2 = ind_j,
            number_pieces = df_pieces_filter$initial_number
          )
        
          dist_ij = dist_ij + distance_pieces * df_pieces_filter$weight
          message(paste0("piece: ", df_pieces_filter$piece,"    colour: ", i_colour,  "   distance result: ", distance_pieces* df_pieces_filter$weight))
        }
      }
      distance_matrix[i,j] = dist_ij
    }
  }
  
  return(distance_matrix)
}

chess_board_example = readr::read_delim(
  file = "data/chess_game/chess_board_example.csv",
  delim = ","
)

if(FALSE){
  View(chess_board_example)
}
      
chess_distance(
  chess_board_example
)
