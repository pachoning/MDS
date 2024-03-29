get_procrustes_parameters <- function(x, target, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (translation) {
    c.target <- scale(target, center = TRUE, scale = FALSE)
    m.target <- attr(c.target ,"scaled:center")

    c.x <- scale(x, center = TRUE, scale = FALSE)
    m.x <- attr(c.x ,"scaled:center")

    matrix_prod <- t(c.target) %*% x
    svd_results <- corpcor::fast.svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)

    translation_vector <- m.target - t(rotation_matrix) %*% m.x

  } else {
    matrix_prod <- t(target) %*% x
    svd_results <- corpcor::fast.svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)
    translation_vector <- matrix(0, n_col, 1)
  }

  return(list(rotation_matrix = rotation_matrix, translation_vector = translation_vector))
}

perform_procrustes <- function(x, target, matrix_to_transform, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (n_row != nrow(target)) {
    stop("x and target do not have same number of rows.\n")
  }

  if (n_col != ncol(target)) {
    stop("x and target do not have same number of columns.\n")
  }

  if (n_col != ncol(matrix_to_transform)) {
    stop("x and matrix_to_transform do not have same number of columns.\n")
  }

  procrustes_parameters <- get_procrustes_parameters(
    x = x,
    target = target,
    translation = translation
  )

  ones_vector <- matrix(
    data = 1,
    nrow = nrow(matrix_to_transform),
    ncol = 1
  )

  translation_matrix <- matrix(
    data = 0,
    nrow = nrow(matrix_to_transform),
    ncol = ncol(matrix_to_transform)
  )
  translation_matrix <- translation_matrix + ones_vector %*% t(procrustes_parameters$translation_vector)
  
  return(matrix_to_transform %*%  procrustes_parameters$rotation_matrix + translation_matrix)
}
