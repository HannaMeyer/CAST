expect_indices_ok <- function(idx, nref) {
  expect_true(is.matrix(idx))
  # allow NA (when k > nref or self excluded) or 1..nref
  vals <- as.vector(idx)
  ok <- is.na(vals) | (vals >= 1 & vals <= nref)
  expect_true(all(ok), info = paste0("indices out of range: ", paste(unique(vals[!ok]), collapse = ", ")))
  # integer-like
  intlike <- sapply(vals[!is.na(vals)], function(x) abs(x - round(x)) < .Machine$double.eps^0.5)
  expect_true(all(intlike), info = "indices are not integer-like")
}

test_that("euclidean: return types and sizes (matrix vs knn)", {
  ref <- as.matrix(iris[1:10, 1:4])
  qry <- as.matrix(iris[11:13, 1:4])
  # full distance matrix
  dm <- knndist(ref, qry, return_distmat = TRUE, dist_fun = "euclidean")
  expect_true(is.matrix(dm))
  expect_type(dm, "double")
  expect_equal(dim(dm), c(nrow(qry), nrow(ref)))
  # k-nearest (k = 2)
  k <- 2
  res <- knndist(ref, qry, k = k, dist_fun = "euclidean", return_distmat = FALSE)
  expect_true(is.matrix(res))
  expect_type(res, "double")
  expect_equal(dim(res), c(nrow(qry), k))
  idx <- attr(res, "indices")
  expect_indices_ok(idx, nrow(ref))
})

test_that("euclidean: single-row query (vector and 1-row data.frame) returns 1 x k", {
  ref <- as.matrix(iris[1:8, 1:4])
  # query as numeric vector
  qvec <- as.numeric(iris[9, 1:4])
  r1 <- knndist(ref, qvec, k = 3, dist_fun = "euclidean")
  expect_true(is.matrix(r1))
  expect_equal(dim(r1), c(1, 3))
  expect_indices_ok(attr(r1, "indices"), nrow(ref))
  # query as 1-row data.frame
  qdf <- iris[9, 1:4, drop = FALSE]
  r2 <- knndist(ref, qdf, k = 2, dist_fun = "euclidean")
  expect_equal(dim(r2), c(1, 2))
  expect_indices_ok(attr(r2, "indices"), nrow(ref))
})

test_that("mahalanobis: basic return types and indices validity", {
  set.seed(101)
  ref <- matrix(rnorm(5 * 3), nrow = 5)
  qry <- matrix(rnorm(2 * 3), nrow = 2)
  res <- knndist(ref, qry, k = 2, dist_fun = "mahalanobis")
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(nrow(qry), 2))
  expect_type(res, "double")
  expect_indices_ok(attr(res, "indices"), nrow(ref))
})

test_that("gower: handles numeric + factor, identical rows -> zero distance to duplicate", {
  ref <- data.frame(
    a = c(0, 1, 0),
    b = factor(c("x", "y", "x")),
    stringsAsFactors = FALSE
  )
  # add an identical duplicate row to guarantee zero non-self distance
  ref <- rbind(ref, ref[1, , drop = FALSE])
  # query NULL: check k = 1 nearest neighbour (should find duplicate with distance 0)
  res <- knndist(ref, query = NULL, k = 1, dist_fun = "gower")
  expect_true(is.matrix(res))
  expect_equal(ncol(res), 1)
  idx <- attr(res, "indices")
  expect_indices_ok(idx, nrow(ref))
  # find positions where duplicate exists: row 1 and row 4 are identical
  # for row 1, nearest neighbour should be row 4 (or vice versa) and distance 0
  expect_true(any(res[, 1] == 0, na.rm = TRUE))
})

test_that("k > n_ref yields NA columns and indices contain NA", {
  ref <- as.matrix(iris[1:5, 1:4])
  qry <- as.matrix(iris[6:7, 1:4])
  k_big <- nrow(ref) + 2
  res <- knndist(ref, qry, k = k_big, dist_fun = "euclidean")
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(nrow(qry), k_big))
  # expect some NA since k > nref
  expect_true(any(is.na(res)))
  idx <- attr(res, "indices")
  expect_true(is.matrix(idx))
  expect_true(any(is.na(idx)))
  # indices that are not NA must be valid
  expect_indices_ok(idx, nrow(ref))
})

test_that("offset affects chosen neighbours", {
  ref <- as.matrix(iris[1:6, 1:4])
  # query is same as reference so self-distances would be NA; use return_distmat=FALSE
  base <- knndist(ref, query = NULL, k = 1, dist_fun = "euclidean", offset = 0)
  off <- knndist(ref, query = NULL, k = 1, dist_fun = "euclidean", offset = 1)
  idx_base <- attr(base, "indices")
  idx_off <- attr(off, "indices")
  # with offset = 1 indices should differ from offset = 0 for at least one row
  expect_true(!all(idx_base == idx_off, na.rm = TRUE))
  expect_indices_ok(idx_off, nrow(ref))
})

test_that("reference single-row and query single-row behave without error", {
  ref <- as.matrix(iris[1, 1:4])
  qry <- as.matrix(iris[2, 1:4])
  res <- knndist(ref, qry, k = 1, dist_fun = "euclidean")
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(1, 1))
  idx <- attr(res, "indices")
  # with only one reference row index should be 1 (or NA if self excluded)
  expect_indices_ok(idx, nrow(ref))
})
