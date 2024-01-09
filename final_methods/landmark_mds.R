# We have rewritten the functions in library lmds to make them 
# comparable (in terms of computing efficiency) to other methods 
landmark_mds <- function (
    x,
    ndim = 3,
    distance_method = c("euclidean", "pearson", "spearman", "cosine", "manhattan"), 
    landmark_method = c("sample"),
    num_landmarks = 500
  ) {
  #dist_2lm <- lmds::select_landmarks(
  dist_2lm <- select_landmarks(
      x = x,
    distance_method = distance_method,
    landmark_method = landmark_method,
    num_landmarks = num_landmarks
  )
  dimred <- base_lmds(dist_2lm = dist_2lm, ndim = ndim, rescale = FALSE)
  rownames(dimred$points) <- rownames(x)
  dimred
}

base_lmds <- function (dist_2lm, ndim = 3, rescale = TRUE, ...){
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
  e <- svd::trlan.eigen(-x_dc/2, neig = ndim)
  ev <- e$d
  evec <- e$u
  ndim1 <- sum(ev > 0)
  if (ndim1 < ndim) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", ndim1, ndim), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
    ndim <- ndim1
  }
  points_inv <- evec/rep(sqrt(ev), each = n)
  #dimred <- -t(t(points_inv) %*% ((dist_2lm - rep(mu_n, each = N))/2))
  #dimred <- (-t(dist_2lm - rep(mu_n, each = N))/2) %*% points_inv
  dimred <- (-t(dist_2lm^2 - rep(mu_n, N))/2) %*% points_inv
  # if (rescale) {
  #   dimred <- dynutils::scale_uniform(dimred)
  # }
  colnames(dimred) <- paste0("comp_", seq_len(ndim))
  dimred <- apply(dimred, MARGIN = 2, FUN = function(y) y - mean(y))
  dimred <- dimred %*% eigen(cov(dimred))$vectors
  mds <- list()
  mds$points <- dimred
  mds$eigen <- ev/n
  return(mds)
}

# We replace lmds::calculate_distance (it is in fact dynutils::calculate_distance) 
# by pracma::distmat
select_landmarks <- function (x, distance_method = c("euclidean", "pearson", "spearman", 
                                                     "cosine", "manhattan"), landmark_method = c("sample"), num_landmarks = 500) 
{
  distance_method <- match.arg(distance_method)
  landmark_method <- match.arg(landmark_method)
  # assert_that(is.matrix(x) || is_sparse(x), is.numeric(num_landmarks), 
  #             length(num_landmarks) == 1, num_landmarks > 2)
  if (num_landmarks > nrow(x)) {
    num_landmarks <- nrow(x)
  }
  if (landmark_method == "sample") {
    ix_lm <- sample.int(nrow(x), num_landmarks)
#    dist_2lm <- as.matrix(calculate_distance(x[ix_lm, , 
#                                               drop = FALSE], x, method = distance_method))
    dist_2lm <- pracma::distmat(X=x[ix_lm, ,drop = FALSE],Y=x)
  }
  #dist_2lm <- zapsmall(dist_2lm)
  attr(dist_2lm, "landmark_ix") <- ix_lm
  dist_2lm
}
