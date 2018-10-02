split_sample <- function(
  x,
  is_random_method = TRUE,
  number_groups = 5,
  variable = NULL,
  return_data_frame = TRUE
){
  
  # Initial parameters ----
  is_data_frame_object = is.data.frame(x)
  is_vector_object = is.vector(x)  
  
  is_numeric = FALSE
  is_date = FALSE
  is_character = FALSE
  is_factor = FALSE
  
  # Convert to a vector
  if(is_data_frame_object == TRUE && is_random_method == FALSE){
    x_transformed = x[[variable]]
    is_numeric = is.numeric(x_transformed)
    is_date = is.Date(x_transformed)
    is_character = is.character(x_transformed)
    is_factor = is.factor(x_transformed)
  } else if(is_data_frame_object == TRUE && is_random_method == TRUE){
    x_transformed = x[[1]]
  } else{
    x_transformed = x
  }
  
  n_obs = length(x_transformed)
  
  # Validations ----
  
  # If it is neither a vector nor a data frame, stop
  if(is_vector_object == FALSE && is_data_frame_object == FALSE){
    stop("x must be either a vector or a data frame")
  }
  
  # If a variable is used, it has to be inside the data frame when a non-radom
  # method is used
  if(
    is_data_frame_object == TRUE  && 
    is_random_method == FALSE   &&
    (
      is.null(variable) |
      any(colnames(x) == variable) == FALSE  
      
    )
  ){
    stop("Please select a valid varible for non-random split method")
  }
  
  # The type of variable has to be a valid one
  if(
    is_random_method == FALSE && 
    any(is_numeric, is_date, is_character, is_factor) == FALSE
  ){
    stop("Invalid type of variable")
  }
  
  # If the number of groups is not well defined, stop
  if(
    (is_numeric == TRUE | is_date == TRUE) && 
    (is.numeric(number_groups) == FALSE | number_groups<1)
  ){
    stop("Define a correct number of number_groups")
  }

  # Random split ----
  if(is_random_method == TRUE){
   
    v_group = rdunif(
      n = n_obs, 
      a = 1,
      b = number_groups 
    )
    
  }


  # Numeric variable ----
  if(is_numeric == TRUE){
    # Getting the quantiles
    labels_quantiles = 1:number_groups
    
    quantiles_to_split = quantcut(
      x = x_transformed,
      q = number_groups,
      labels = labels_quantiles
    )
    
    v_group = quantiles_to_split
  }
  
  # Date variable  ----
  if(is_date == TRUE){
    # To get the quantiles, the date is transformed to numeric
    date_to_numeric = as.numeric(x_transformed) 
    
    # Getting the quantiles
    labels_quantiles = 1:number_groups
    
    quantiles_to_split = quantcut(
      x = date_to_numeric,
      q = number_groups,
      labels = labels_quantiles
    )
      
    v_group = quantiles_to_split
      
  }
  
  
  # String variable ----
  if( is_character == TRUE ){
    x_factor = as.factor(x_transformed)
    levels(x_factor) = 1:length(levels(x_factor))
    x_factor = as.numeric(x_factor)
    
    v_group = x_factor
  }
  
  
  # Factor variable ----
  if(is_factor == TRUE){
    x_factor = as.factor(x_transformed)
    levels(x_factor) = 1:length(levels(x_factor))
    x_factor = as.numeric(x_factor)
    
    v_group = x_factor
  }

  
  # output ----
  if( return_data_frame == TRUE){
    return(
      cbind(
        x,
        group_member = v_group
      )
    )
  } else{
    return(
      split(
        x = x,
        f = v_group
      )
    )
  }
  
}
