#' Compute k-nearest neighbor distances
#'
#' Compute distances between rows of a reference dataset and (optionally) a query
#' dataset and return the k-nearest neighbor distances (with neighbor indices).
#'
#' @param reference A numeric matrix, numeric vector, or data.frame of reference
#'   observations (rows = observations, columns = features).
#' @param query A numeric matrix, numeric vector, data.frame of query
#'   observations, or NULL. If NULL (default) distances are computed among rows
#'   of \code{reference}.
#' @param k Integer; number of neighbours to return (default 1).
#' @param dist_fun Character; distance metric. One of \code{"euclidean"},
#'   \code{"mahalanobis"} or \code{"gower"} (default \code{"euclidean"}).
#'   - \code{"gower"}: numeric columns of \code{reference} are scaled to [0,1]
#'     using the reference min/max; categorical columns (character/factor)
#'     are converted to integer codes before distance calculation.
#'   - \code{"mahalanobis"}: uses a pseudo-inverse covariance (MASS::ginv) and
#'     transforms data so that Euclidean distance in transformed space equals
#'     Mahalanobis distance in the original space.
#' @param offset Integer offset for neighbor ranking (default 0). The returned
#'   neighbors are selected from the sorted distance vector using the range
#'   \code{c(max(1, 1 + offset), max(k, offset + k))}. This can be used to skip the
#'   immediate nearest neighbour(s) when desired, e.g. to exclude the observation itself when
#'   \code{query = NULL}.
#' @param return_distmat Logical; if \code{TRUE} returns the full distance
#'   matrix (matrix with rows = query, cols = reference). If \code{FALSE}
#'   (default) returns the k-nearest distances matrix with an \code{"indices"}
#'   attribute (see Value).
#'
#' @return If \code{return_distmat = TRUE} a numeric distance matrix is
#'   returned. If \code{return_distmat = FALSE} a numeric matrix of nearest
#'   neighbor distances is returned (rows correspond to query observations,
#'   columns correspond to the requested neighbour ranks). The returned object
#'   has an attribute \code{"indices"} containing an integer matrix of the
#'   corresponding neighbour row indices in \code{reference}.
#'
#' @details
#' - When \code{query = NULL} self-distances are set to \code{NA} so that the
#'   nearest neighbour search excludes the observation itself.
#' - Categorical variables are converted to integer factor codes only when
#'   \code{dist_fun = "gower"} (via an internal helper).
#' - For \code{"mahalanobis"} a numerical linear-algebra fallback using eigen
#'   decomposition is used if Cholesky fails.
#'
#' @examples
#' # nearest 2 neighbours within iris (euclidean)
#' res <- CAST:::knndist(iris[,1:4], k = 2)
#' str(res)
#' # access neighbour indices
#' attr(res, "indices")[1, ]
#'
#' # distances from first 5 rows of iris to the full reference
#' CAST:::knndist(iris[,1:4], iris[1:5,1:4], k = 3)
#'
#' # return full distance matrix
#' dm <- CAST:::knndist(iris[,1:4], return_distmat = TRUE)
#' dim(dm)
#'
#' @seealso \code{\link[philentropy]{distance}},
#'   \code{\link[philentropy]{dist_many_many}}
#' @keywords internal
knndist <- function(
  reference,
  query = NULL,
  k = 1,
  dist_fun = c("euclidean", "mahalanobis", "gower"),
  offset = 0,
  return_distmat = FALSE
) {
  dist_fun <- match.arg(dist_fun)

  if (dist_fun == "gower") {
    gower_res <- .preprocess_gower(reference, query)
    reference <- gower_res$reference
    query <- gower_res$query
  }

  # normalize input to matrices (handle numeric vectors and data.frames)
  if (inherits(reference, "numeric")) {
    reference <- matrix(reference, nrow = 1)
  }
  if (inherits(reference, "data.frame")) {
    reference <- as.matrix(reference)
  }
  if (inherits(query, "numeric")) {
    query <- matrix(query, nrow = 1)
  }
  if (inherits(query, "data.frame")) {
    query <- as.matrix(query)
  }

  if (dist_fun == "mahalanobis") {
    maha_res <- .preprocess_maha(reference, query)
    reference <- maha_res$reference
    query <- maha_res$query
    dist_fun <- "euclidean"
  }

  dist_mat <- .dist_mat(reference, query, dist_fun)

  if (return_distmat) {
    return(dist_mat)
  }

  knn_dists <- .knn_dist(dist_mat, k, offset)
  return(knn_dists)

}

.dist_mat <- function(reference, query = NULL, dist_fun) {
  diag_to_na <- if(is.null(query)) TRUE else FALSE
  if (is.null(query)) query <- reference
  dist_mat <- philentropy::dist_many_many(dists1 = query, dists2 = reference, method = dist_fun)
  if (diag_to_na) diag(dist_mat) <- NA_real_
  return(dist_mat)
}

.knn_dist <- function(dist_mat, k, offset) {
  range <- c(max(1, 1 + offset), max(k, offset + k))
  nc <- diff(range) + 1L
  nr <- nrow(dist_mat)
  # preallocate
  knn_dists <- matrix(NA_real_, nrow = nr, ncol = nc)
  knn_indices <- matrix(NA_integer_, nrow = nr, ncol = nc)
  # loop by row (order with NA last so NA self-distances are excluded)
  for (i in seq_len(nr)) {
    row <- dist_mat[i, ]
    o <- order(row, seq_along(row), na.last = TRUE)
    idx <- o[range[1]:range[2]]
    knn_indices[i, ] <- as.integer(idx)
    knn_dists[i, ] <- row[idx]
  }
  attr(knn_dists, "indices") <- knn_indices
  return(knn_dists)
}

.preprocess_maha <- function(reference, query = NULL) {
  # For Mahalanobis distance, we need to compute the inverse covariance matrix
  # we then transform that reference and query to calculate the L2 distance in
  # the transformed space, which is equivalent to the Mahalanobis distance in the original space.
  S_inv <- MASS::ginv(stats::cov(reference))
  chol_ok <- try(R <- chol(S_inv), silent = TRUE)
  if (!inherits(chol_ok, "try-error")) {
    A <- t(R)
  } else {
    eig <- eigen(S_inv, symmetric = TRUE)
    vals <- pmax(eig$values, 0) # guard against tiny negative values
    A <- eig$vectors %*% diag(sqrt(vals)) %*% t(eig$vectors)
  }
  reference <- reference %*% A
  if (!is.null(query)) {
    query <- query %*% A
  }
  return(list(reference = reference, query = query))
}

.preprocess_gower <- function(reference, query = NULL) {
  # For Gower distance, we need to scale numeric columns to [0,1] and convert
  # categorical columns to integer codes. We use the reference for scaling and
  # encoding to ensure consistency between reference and query.
  is_cat <- sapply(reference, function(col) is.character(col) || is.factor(col))
  if (any(is_cat)) {
    reference[is_cat] <- lapply(reference[is_cat], function(col) as.integer(as.factor(col)))
    if (!is.null(query)) {
      query[is_cat] <- lapply(query[is_cat], function(col) as.integer(as.factor(col)))
    }
  }
  is_num <- !is_cat
  if (any(is_num)) {
    mins <- apply(reference[is_num], 2, min)
    maxs <- apply(reference[is_num], 2, max)
    ranges <- maxs - mins
    # avoid division by zero for constant columns
    ranges[ranges == 0] <- 1
    reference[is_num] <- sweep(sweep(reference[is_num], 2, mins), 2, ranges, FUN = "/")
    if (!is.null(query)) {
      query[is_num] <- sweep(sweep(query[is_num], 2, mins), 2, ranges, FUN = "/")
    }
  }
  return(list(reference = reference, query = query))
}
