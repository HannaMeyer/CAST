#' Compute pairwise distances between reference and query
#'
#' Compute a pairwise distance matrix between rows of reference and (optionally)
#' query using the selected metric. This function performs necessary
#' preprocessing for some metrics (Gower scaling, Mahalanobis transformation)
#' and then calls philentropy::dist_many_many to compute distances.
#'
#' @param reference A numeric matrix, numeric vector, or data.frame of reference
#'   observations (rows = observations, columns = features). For Gower and
#'   Mahalanobis preprocessing the reference is used to compute scaling/transform
#'   parameters applied to query (if provided).
#' @param query A numeric matrix, numeric vector, data.frame of query
#'   observations, or NULL. If NULL (default) distances are computed among rows
#'   of reference and the self-distance is set to NA so an observation does not
#'   return itself as a nearest neighbour.
#' @param dist_fun Character; distance metric to use. One of
#'   "euclidean", "mahalanobis" or "gower" (default "euclidean"). See
#'   Details for preprocessing performed for non-euclidean metrics.
#'
#' @return A numeric distance matrix with rows = query and cols = reference. If
#'   query is NULL the distance matrix is computed among reference rows and the
#'   diagonal is set to NA so self-distances are excluded.
#'
#' @details
#' - For "gower" numeric reference columns are scaled to [0,1] using reference
#'   min/max and categorical columns are converted to integer factor codes. The
#'   same scaling/encoding is applied to query (when provided) using the
#'   reference-derived parameters.
#' - For "mahalanobis" the reference covariance is pseudo-inverted (MASS::ginv)
#'   and a linear transform is applied so that Euclidean distance in transformed
#'   space equals Mahalanobis distance in the original space. A numerical
#'   eigenvalue fallback is used if Cholesky fails.
#'
#' @seealso \code{\link{.knndist}}, \code{\link[philentropy]{dist_many_many}}
#' @keywords internal
.distance <- function(
  reference,
  query = NULL,
  dist_fun = c("euclidean", "mahalanobis", "gower")
) {
  dist_fun <- match.arg(dist_fun)

  preprocess <- switch(
    dist_fun,
    "mahalanobis" = .preprocess_maha,
    "gower" = .preprocess_gower,
    .preprocess_default # for "euclidean" and other metrics
  )

  pre_res <- preprocess(reference, query, dist_fun)
  dist_mat <- do.call(.dist_mat, pre_res)
  return(dist_mat)
}

#' Compute k-nearest neighbor distances
#'
#' Compute distances between rows of a reference dataset and (optionally) a query
#' dataset and return the k-nearest neighbor distances (with neighbour indices).
#'
#' This function delegates pairwise distance computation to \code{\link{.distance}} and then
#' selects the k nearest neighbours per query row.
#'
#' @inheritParams .distance
#' @param k Integer; number of neighbours to return (default 1).
#' @param offset Integer offset for neighbor ranking (default 0). The returned
#'   neighbours are selected from the sorted distance vector using the range
#'   c(max(1, 1 + offset), max(k, offset + k)). This can be used to skip the
#'   immediate nearest neighbour(s) when desired, e.g. to exclude the
#'   observation itself when query = NULL.
#'
#' @return A numeric matrix of nearest neighbour distances (rows correspond to
#'   query observations, columns correspond to neighbour ranks). The returned
#'   object has an attribute "indices" containing an integer matrix of the
#'   corresponding neighbour row indices in reference.
#'
#' @details
#' knndist delegates distance computation to \code{\link{.distance}} and then performs the
#' neighbour ranking. For behaviour of the individual distance metrics and any
#' preprocessing (e.g. scaling for Gower or Mahalanobis transformation) see
#' \code{\link{.distance}}.
#'
#' @seealso \code{\link{.distance}}, \code{\link[philentropy]{dist_many_many}}
#' @keywords internal
.knndist <- function(
  reference,
  query = NULL,
  k = 1,
  dist_fun = c("euclidean", "mahalanobis", "gower"),
  offset = 0
) {
  dist_mat <- .distance(reference, query, dist_fun)
  knn_dists <- .knn_dist(dist_mat, k, offset)
  return(knn_dists)
}

################################################################################
# Core distance matrix computation and k-nearest neighbor selection
################################################################################

.dist_mat <- function(reference, query = NULL, dist_fun) {
  diag_to_na <- if (is.null(query)) TRUE else FALSE
  if (is.null(query)) {
    query <- reference
  }
  dist_mat <- philentropy::dist_many_many(
    dists1 = query,
    dists2 = reference,
    method = dist_fun
  )
  if (diag_to_na) {
    diag(dist_mat) <- NA_real_
  }
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

################################################################################
# Data normalization and preprocessing
################################################################################

.normalize_to_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  # normalize input to matrices (handle numeric vectors and data.frames)
  if (inherits(x, "numeric")) {
    x <- matrix(x, nrow = 1)
  }
  if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  }
  return(x)
}


# preprocessing is expected to return a list with elements "reference", "query", and "dist_fun"
# and handle matrix normalization (e.g. from data.frame to matrix) as needed
.preprocess_default <- function(ref, qry, fun = "euclidean") {
  ref <- .normalize_to_matrix(ref)
  query <- .normalize_to_matrix(qry)
  return(list(reference = ref, query = query, dist_fun = fun))
}

# For Mahalanobis distance, we need to compute the inverse covariance matrix
# we then transform that reference and query to calculate the L2 distance in
# the transformed space, which is equivalent to the Mahalanobis distance in the original space.
.preprocess_maha <- function(reference, query = NULL, dist_fun = "mahalanobis") {
  reference <- .normalize_to_matrix(reference)
  query <- .normalize_to_matrix(query)

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
  reference <- .normalize_to_matrix(reference)
  query <- .normalize_to_matrix(query)
  return(list(reference = reference, query = query, dist_fun = "euclidean"))
}

  # For Gower distance, we need to scale numeric columns to [0,1] and convert
# categorical columns to integer codes. We use the reference for scaling and
# encoding to ensure consistency between reference and query.

.preprocess_gower <- function(reference, query = NULL, dist_fun = "gower") {
  reference <- .normalize_to_df(reference)
  query <- .normalize_to_df(query)
  has_query <- !is.null(query)

  # drop unknown levels and set unused levels to NA in reference and query for each categorical variable
  catvars <- .get_categorical_variables(reference)

  if (length(catvars) > 0) {
    for (catvar in catvars) {
      res <- .drop_unknown_levels(reference, query, catvar)
      # TODO: rename outputs of drop_unknown_levels to be more intuitive (e.g. reference and query instead of train and newdata)
      reference <- res$reference
      reference[[catvar]] <- as.integer(as.factor(reference[[catvar]]))
      query <- res$query
      if (has_query) {
        query[[catvar]] <- as.integer(as.factor(query[[catvar]]))
      }
    }
  }

  # numeric variable scaling to [0,1] using reference min/max
  numvars <- setdiff(names(reference), catvars)

  if (length(numvars) > 0) {
    mins <- apply(reference[numvars], 2, min)
    maxs <- apply(reference[numvars], 2, max)

    reference[numvars] <- .minmax_normalize(reference, numvars, maxs, mins)
    if (has_query) {
      query[numvars] <- .minmax_normalize(query, numvars, maxs, mins)
    }
  }

  # convert back to matrices for distance computation
  reference <- .normalize_to_matrix(reference)
  query <- .normalize_to_matrix(query)

  return(list(reference = reference, query = query, dist_fun = "gower"))
}

# matrix -> df; vector -> df with one row; df -> df (no change)
.normalize_to_df <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  } else if (inherits(x, "numeric")) {
    x <- as.data.frame(t(x))
  }
  return(x)
}

# helper for Gower scaling numeric columns to [0,1] using reference min/max
.minmax_normalize <- function(df, vars, maxs, mins) {
  ranges <- maxs - mins
  # avoid division by zero for constant columns
  ranges[ranges == 0] <- 1
  return(sweep(sweep(df[vars], 2, mins), 2, ranges, FUN = "/"))
}