#procrustes_based_distance <- function(x, target, list_dilation_factor, list_rotation_matrix, list_translation_vector) {
#  
#  ones_vector = matrix(data = 1, nrow = nrow(target), ncol = 1)
#  translation_matrix = matrix(data = 0, nrow = nrow(target), ncol = ncol(target))
#  translation_matrix = lapply(list_translation_vector, 
#                              function(tv, tm, ov) tm + ov %*% t(tv),
#                              tm = translation_matrix, ov = ones_vector)
#  
#  proc_trans = mapply(function(dil_fac, rot_m, trans_m, x) dil_fac * x %*% rot_m + trans_m,
#                      dil_fac = list_dilation_factor, 
#                      rot_m = list_rotation_matrix, 
#                      trans_m = translation_matrix,
#                      MoreArgs = list(x = x),
#                      SIMPLIFY = FALSE)
#  
#  frob_norm = sapply(proc_trans, function(x, y) norm(x -y, type = "F"), y = target)
#  
#  best_norm = which.min(frob_norm)
#  
#  return(list(dilation_factor = list_dilation_factor[[best_norm]],
#              rotation_matrix = list_rotation_matrix[[best_norm]],
#              translation_vector = list_translation_vector[[best_norm]]))
#
#}



procrustes_based_distance <- function(list_dilation_factor, list_rotation_matrix, list_translation_vector) {

  total_elements = length(list_rotation_matrix)
  dist_matrix = matrix(data = NA, nrow = total_elements, total_elements)

  for (i in 1:total_elements) {
    current_matrix = list_rotation_matrix[[i]]
    for (j in 1:total_elements) {
      other_matrix = list_rotation_matrix[[j]]
      diff_matrix = current_matrix - other_matrix
      dist_matrix[i, j] = norm(diff_matrix, type = "F")
    }
  }

  if (total_elements > 1) {
    dist_filtered = sapply(1:total_elements, function(ind, matrix) matrix[ind, -ind], matrix = dist_matrix)
    median_distance = apply(X = dist_matrix, MARGIN = 1, FUN = median, na.rm = TRUE)
    arg_min_median_distance = which.min(median_distance)
  } else {
    arg_min_median_distance = 1
  }

  return(list(dilation_factor = list_dilation_factor[[arg_min_median_distance]],
              rotation_matrix = list_rotation_matrix[[arg_min_median_distance]],
              translation_vector = list_translation_vector[[arg_min_median_distance]]))

}


get_procrustes_parameters <- function(x, target, translation = FALSE, dilation = FALSE) {

  n_row = nrow(x)
  n_col = ncol(x)

  diag_matrix = diag(n_row)
  if (translation) {
    diag_matrix = diag(n_row) - 1/n_row * matrix(1, n_row, n_row)
  }

  matrix_prod = t(target) %*% diag_matrix %*% x
  svd_results = svd(matrix_prod)
  rotation_matrix = svd_results$v %*% t(svd_results$u)

  dilation_factor = 1
  if (dilation) {
    mat1 = t(target) %*% diag_matrix %*% x %*% rotation_matrix
    mat2 = t(x) %*% diag_matrix %*% x
    num = 0
    denom = 0

    for (i in 1:n_col) {
      num = num + mat1[i, i]
      denom = denom + mat2[i, i]
    }

    dilation_factor = num/denom
  }

  translation_vector = matrix(0, n_col, 1)
  if (translation) {
    translation_vector = 1/n_row * t(target - dilation_factor * x %*% rotation_matrix) %*% matrix(1, n_row, 1)
  }

  return(list(dilation_factor = dilation_factor, rotation_matrix = rotation_matrix, translation_vector = translation_vector))
}


perform_procrustes <- function(x, target, matrix_to_transform, translation = FALSE, dilation = FALSE, 
                               largest_matrix_efficient_procrustes = 10000) {

  n_row = nrow(x)
  n_col = ncol(x)

  if (n_row != nrow(target)) {
    stop("x and target do not have same number of rows.\n")
  }

  if (n_col != ncol(target)) {
    stop("x and target do not have same number of columns.\n")
  }

  if (n_col != ncol(matrix_to_transform)) {
    stop("x and matrix_to_transform do not have same number of columns.\n")
  }

  p = 1
  if (nrow(x)>largest_matrix_efficient_procrustes) {
    p = ceiling(nrow(x)/largest_matrix_efficient_procrustes) + 2
  }

  indexes_group = sample(x = p, size = nrow(x), replace = TRUE)
  indexes_group = lapply(1:p, function(x, y) which(x==y), y = indexes_group)

  x_filtered = lapply(indexes_group, function(idx, matrix) matrix[idx, ,drop = FALSE], matrix = x)
  target_filtered = lapply(indexes_group, function(idx, matrix) matrix[idx, ,drop = FALSE], matrix = target)

  multiple_procrustes_results = mapply(get_procrustes_parameters, 
                                       x = x_filtered, 
                                       target = target_filtered,
                                       MoreArgs = list(translation = translation, dilation = dilation), 
                                       SIMPLIFY = FALSE)

  list_dilation_factor = lapply(multiple_procrustes_results, function(x) x$dilation_factor)  
  list_rotation_matrix = lapply(multiple_procrustes_results, function(x) x$rotation_matrix)
  list_translation_vector = lapply(multiple_procrustes_results, function(x) x$translation)

  procrustes_parameters = procrustes_based_distance(list_dilation_factor = list_dilation_factor,
                                                    list_rotation_matrix = list_rotation_matrix,
                                                    list_translation_vector = list_translation_vector)

  ones_vector = matrix(data = 1, nrow=nrow(matrix_to_transform), ncol = 1)
  translation_matrix = matrix(data = 0, nrow = nrow(matrix_to_transform), ncol = ncol(matrix_to_transform))
  translation_matrix = translation_matrix + ones_vector %*% t(procrustes_parameters$translation_vector)

  return(procrustes_parameters$dilation_factor * matrix_to_transform %*%  procrustes_parameters$rotation_matrix + translation_matrix)
}
