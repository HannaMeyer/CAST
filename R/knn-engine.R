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
#' res <- knndist(iris[,1:4], k = 2)
#' str(res)
#' # access neighbour indices
#' attr(res, "indices")[1, ]
#'
#' # distances from first 5 rows of iris to the full reference
#' knndist(iris[,1:4], iris[1:5,1:4], k = 3)
#'
#' # return full distance matrix
#' dm <- knndist(iris[,1:4], return_distmat = TRUE)
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
  return_distmat = FALSE) {
  
  dist_fun <- match.arg(dist_fun)

  if (dist_fun == "gower") {
    # requires scaling of the numerical variables, while categorical variables are not scaled.
    num_vars <- vapply(reference, is.numeric, logical(1))
    ref_num <- reference[, num_vars, drop = FALSE]
    mins <- sapply(ref_num, min, na.rm = TRUE)
    maxs <- sapply(ref_num, max, na.rm = TRUE)
    range <- maxs - mins

    reference[, num_vars] <- as.data.frame(
        sweep(sweep(reference[, num_vars, drop = FALSE], 2, mins, "-"), 2, range, "/"),
        stringsAsFactors = FALSE
    )
    
    if (!is.null(query)) {
      query[, num_vars] <- as.data.frame(
        sweep(sweep(query[, num_vars, drop = FALSE], 2, mins, "-"), 2, range, "/"),
        stringsAsFactors = FALSE
      )
    }
    query <- .numeric_fct(query)
    reference <- .numeric_fct(reference)
  }
  
  if (inherits(query, "numeric")) {
    query <- matrix(query, nrow = 1)
  }
  if (inherits(query, "data.frame")) {
    query <- as.matrix(query)
  }
  if (inherits(reference, "numeric")) {
    reference <- matrix(reference, nrow = 1)
  }
  if (inherits(reference, "data.frame")) {
    reference <- as.matrix(reference)
  }

  if (dist_fun == "mahalanobis") {
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
    if (!is.null(query)) query <- query %*% A
    dist_fun = "euclidean"
  }


  # calculate the distance matrix
  if (is.null(query)) {
    dists <- suppressMessages(philentropy::distance(reference, method = dist_fun))
    if (length(dists) == 1) return(dists)
    diag(dists) <- NA # Exclude self-distance
    } else {
    dists <- suppressMessages(philentropy::dist_many_many(query, reference, method = dist_fun))
  }

  if (return_distmat) return(dists)

  # compute range and sizes
  range <- c(max(1, 1 + offset), max(k, offset + k))
  nc <- diff(range) + 1L; nr <- nrow(dists)

  # preallocate
  knn_dists <- matrix(NA_real_, nrow = nr, ncol = nc)
  knn_indices <- matrix(NA_integer_, nrow = nr, ncol = nc)

  # loop by row (order with NA last so NA self-distances are excluded)
  for (i in seq_len(nr)) {
    row <- dists[i, ]
    o <- order(row, na.last = TRUE)
    idx <- o[range[1]:range[2]]
    knn_indices[i, ] <- as.integer(idx)
    knn_dists[i, ] <- row[idx]
  }

  attr(knn_dists, "indices") <- knn_indices
  knn_dists

}  


#' Convert character/factor columns to integer codes (internal)
#'
#' Helper to convert factor/character columns of a data.frame to integer codes
#' (as.integer(as.factor(.))). Returns the input unchanged if there are no
#' character/factor columns. Used when computing Gower distances.
#'
#' @param x A data.frame (or \code{NULL}).
#' @return The input data.frame with any factor/character columns replaced by
#'   integer codes, or \code{NULL} if \code{x} is \code{NULL}.
#' @keywords internal
#' @noRd
.numeric_fct <- function(x) {
  if (is.null(x)) return(NULL)
  catVars <- names(x)[vapply(x, function(z) inherits(z, c("factor", "character")), logical(1))]
  if (length(catVars) == 0) return(x)
  # Convert categorical variables to factor then to numeric
  for (var in catVars) {
    x[[var]] <- as.integer(as.factor(x[[var]]))
  }
  return(x)
}