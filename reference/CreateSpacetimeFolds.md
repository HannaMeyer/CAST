# Create Space-time Folds

Create spatial, temporal or spatio-temporal Folds for cross validation
based on pre-defined groups

## Usage

``` r
CreateSpacetimeFolds(
  x,
  spacevar = NA,
  timevar = NA,
  k = 10,
  class = NA,
  seed = sample(1:1000, 1)
)
```

## Arguments

- x:

  data.frame containing spatio-temporal data

- spacevar:

  Character indicating which column of x identifies the spatial units
  (e.g. ID of weather stations)

- timevar:

  Character indicating which column of x identifies the temporal units
  (e.g. the day of the year)

- k:

  numeric. Number of folds. If spacevar or timevar is NA and a leave one
  location out or leave one time step out cv should be performed, set k
  to the number of unique spatial or temporal units.

- class:

  Character indicating which column of x identifies a class unit (e.g.
  land cover)

- seed:

  numeric. See ?seed

## Value

A list that contains a list for model training and a list for model
validation that can directly be used as "index" and "indexOut" in
caret's trainControl function. "cluster" gives us the information to
which validation fold a sample belongs.

## Details

The function creates train and test sets by taking (spatial and/or
temporal) groups into account. In contrast to
[`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md), it
requires that the groups are already defined (e.g. spatial clusters or
blocks or temporal units). Using "class" is helpful in the case that
data are clustered in space and are categorical. E.g This is the case
for land cover classifications when training data come as training
polygons. In this case the data should be split in a way that entire
polygons are held back (spacevar="polygonID") but at the same time the
distribution of classes should be similar in each fold (class="LUC").

## Note

Standard k-fold cross-validation can lead to considerable
misinterpretation in spatial-temporal modelling tasks. This function can
be used to prepare a Leave-Location-Out, Leave-Time-Out or
Leave-Location-and-Time-Out cross-validation as target-oriented
validation strategies for spatial-temporal prediction tasks. See Meyer
et al. (2018) for further information. CreateSpaceTimeFolds is just a
very simple approach and the suitability depends on the choice of the
groups. You may check the suitability with
[`geodist`](https://hannameyer.github.io/CAST/reference/geodist.md).
Consider [`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md)
or [`knndm`](https://hannameyer.github.io/CAST/reference/knndm.md) as
alternatives or other approaches such as Spatial Blocks. For spatial
visualization of fold affiliation see examples.

## References

Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018):
Improving performance of spatio-temporal machine learning models using
forward feature selection and target-oriented validation. Environmental
Modelling & Software 101: 1-9.

## See also

[`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html),[`ffs`](https://hannameyer.github.io/CAST/reference/ffs.md),
[`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md),
[`geodist`](https://hannameyer.github.io/CAST/reference/geodist.md)

## Author

Hanna Meyer

## Examples

``` r
if (FALSE) { # \dontrun{
data(cookfarm)
### Prepare for 10-fold Leave-Location-and-Time-Out cross validation
indices <- CreateSpacetimeFolds(cookfarm,"SOURCEID","Date")
str(indices)
### Prepare for 10-fold Leave-Location-Out cross validation
indices <- CreateSpacetimeFolds(cookfarm,spacevar="SOURCEID")
str(indices)
### Prepare for leave-One-Location-Out cross validation
indices <- CreateSpacetimeFolds(cookfarm,spacevar="SOURCEID",
    k=length(unique(cookfarm$SOURCEID)))
str(indices)

### example from splotopen and visualization
data(splotdata)
indices <- CreateSpacetimeFolds(splotdata,spacevar="Country")
ggplot() +
geom_sf(data = splotdata, aes(col = factor(indices$cluster)))
## is this representative?
data(splotdata)
studyArea <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf")
dist <- geodist(splotdata, studyArea,cvfolds=indices$cluster)
plot(dist)+ scale_x_log10(labels=round)

} # }
```
