# Compute pairwise distances between reference and query

Compute a pairwise distance matrix between rows of reference and
(optionally) query using the selected metric. This function performs
necessary preprocessing for some metrics (Gower scaling, Mahalanobis
transformation) and then calls philentropy::dist_many_many to compute
distances.

## Usage

``` r
.distance(
  reference,
  query = NULL,
  dist_fun = c("euclidean", "mahalanobis", "gower")
)
```

## Arguments

- reference:

  A numeric matrix, numeric vector, or data.frame of reference
  observations (rows = observations, columns = features). For Gower and
  Mahalanobis preprocessing the reference is used to compute
  scaling/transform parameters applied to query (if provided).

- query:

  A numeric matrix, numeric vector, data.frame of query observations, or
  NULL. If NULL (default) distances are computed among rows of reference
  and the self-distance is set to NA so an observation does not return
  itself as a nearest neighbour.

- dist_fun:

  Character; distance metric to use. One of "euclidean", "mahalanobis"
  or "gower" (default "euclidean"). See Details for preprocessing
  performed for non-euclidean metrics.

## Value

A numeric distance matrix with rows = query and cols = reference. If
query is NULL the distance matrix is computed among reference rows and
the diagonal is set to NA so self-distances are excluded.

## Details

\- For "gower" numeric reference columns are scaled to \[0,1\] using
reference min/max and categorical columns are converted to integer
factor codes. The same scaling/encoding is applied to query (when
provided) using the reference-derived parameters. - For "mahalanobis"
the reference covariance is pseudo-inverted (MASS::ginv) and a linear
transform is applied so that Euclidean distance in transformed space
equals Mahalanobis distance in the original space. A numerical
eigenvalue fallback is used if Cholesky fails.

## See also

[`.knndist`](https://hannameyer.github.io/CAST/reference/dot-knndist.md),
[`dist_many_many`](https://drostlab.github.io/philentropy/reference/dist_many_many.html)
