my_lmds <- function (x, ndim = 3, distance_method = c("euclidean", "pearson", 
                                           "spearman", "cosine", "manhattan"), 
                     landmark_method = c("sample"), 
          num_landmarks = 500) 
{
  require(lmds)
  dist_2lm <- lmds::select_landmarks(x = x, distance_method = distance_method, 
                               landmark_method = landmark_method, num_landmarks = num_landmarks)
  dimred <- my_cmdscale_landmarks(dist_2lm = dist_2lm, ndim = ndim, 
                               rescale = TRUE)
  rownames(dimred) <- rownames(x)
  dimred
}


my_cmdscale_landmarks <- function (dist_2lm, ndim = 3, rescale = TRUE, ...) 
{
  # assert_that(!is.null(attr(dist_2lm, "landmark_ix")) || (!is.null(rownames(dist_2lm)) && 
  #                                                           !is.null(colnames(dist_2lm))), is.numeric(ndim), length(ndim) == 
  #               1, ndim >= 1, is.logical(rescale), length(rescale) == 
  #               1, !is.na(rescale))
  ix_lm <- attr(dist_2lm, "landmark_ix")
  if (is.null(ix_lm)) {
    ix_lm <- match(rownames(dist_2lm), colnames(dist_2lm))
  }
  x <- dist_2lm[, ix_lm, drop = FALSE]^2
  n <- as.integer(nrow(x))
  N <- as.integer(ncol(dist_2lm))
  mu_n <- rowMeans(x)
  mu <- mean(x)
  x_c <- sweep(x, 1, mu_n, "-")
  x_dc <- sweep(x_c, 2, mu_n, "-") + mu
  # if (ndim > 0.5 * min(nrow(x_dc), ncol(x_dc))) {
  #   e <- eigen(-x_dc/2, symmetric = TRUE)
  #   e$values <- e$values[seq_len(ndim)]
  #   e$vectors <- e$vectors[, seq_len(ndim), drop = FALSE]
  # }
  # else {
  #   e <- irlba::partial_eigen(-x_dc/2, symmetric = TRUE, 
  #                             n = ndim, ...)
  # }
  # ev <- e$values
  # evec <- e$vectors
  e <- svd::trlan.eigen(-x_dc/2,neig = ndim)
  ev <- e$d
  evec <- e$u
  ndim1 <- sum(ev > 0)
  if (ndim1 < ndim) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", 
                     ndim1, ndim), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
    ndim <- ndim1
  }
  points_inv <- evec/rep(sqrt(ev), each = n)
  #dimred <- -t(t(points_inv) %*% ((dist_2lm - rep(mu_n, each = N))/2))
  dimred <- (-t(dist_2lm - rep(mu_n, each = N))/2) %*% points_inv
  if (rescale) {
    dimred <- dynutils::scale_uniform(dimred)
  }
  colnames(dimred) <- paste0("comp_", seq_len(ndim))
  dimred
}
