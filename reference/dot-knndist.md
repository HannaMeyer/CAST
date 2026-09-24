# Compute k-nearest neighbor distances

Compute distances between rows of a reference dataset and (optionally) a
query dataset and return the k-nearest neighbor distances (with
neighbour indices).

## Usage

``` r
.knndist(
  reference,
  query = NULL,
  k = 1,
  dist_fun = c("euclidean", "mahalanobis", "gower"),
  offset = 0
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

- k:

  Integer; number of neighbours to return (default 1).

- dist_fun:

  Character; distance metric to use. One of "euclidean", "mahalanobis"
  or "gower" (default "euclidean"). See Details for preprocessing
  performed for non-euclidean metrics.

- offset:

  Integer offset for neighbor ranking (default 0). The returned
  neighbours are selected from the sorted distance vector using the
  range c(max(1, 1 + offset), max(k, offset + k)). This can be used to
  skip the immediate nearest neighbour(s) when desired, e.g. to exclude
  the observation itself when query = NULL.

## Value

A numeric matrix of nearest neighbour distances (rows correspond to
query observations, columns correspond to neighbour ranks). The returned
object has an attribute "indices" containing an integer matrix of the
corresponding neighbour row indices in reference.

## Details

This function delegates pairwise distance computation to
[`.distance`](https://hannameyer.github.io/CAST/reference/dot-distance.md)
and then selects the k nearest neighbours per query row.

knndist delegates distance computation to
[`.distance`](https://hannameyer.github.io/CAST/reference/dot-distance.md)
and then performs the neighbour ranking. For behaviour of the individual
distance metrics and any preprocessing (e.g. scaling for Gower or
Mahalanobis transformation) see
[`.distance`](https://hannameyer.github.io/CAST/reference/dot-distance.md).

## See also

[`.distance`](https://hannameyer.github.io/CAST/reference/dot-distance.md),
[`dist_many_many`](https://drostlab.github.io/philentropy/reference/dist_many_many.html)
