# Calculate euclidean nearest neighbor distances in geographic space or feature space

Calculates nearest neighbor distances in geographic space or feature
space between training data as well as between prediction locations and
training data. Optional, the nearest neighbor distances between test
data and training data or between different CV folds is computed.

## Usage

``` r
geodist(
  x,
  modeldomain = NULL,
  dist_space = "geographical",
  CVtest = NULL,
  CVtrain = NULL,
  testdata = NULL,
  preddata = NULL,
  samplesize = 2000,
  sampling = "regular",
  variables = NULL,
  time_var = NULL,
  time_unit = "auto",
  algorithm = "brute",
  dist_fun = "euclidean",
  scale_vars = TRUE,
  cvtrain = NULL,
  cvfolds = NULL,
  type = NULL,
  timevar = NULL
)
```

## Arguments

- x:

  object of class sf, training data locations

- modeldomain:

  SpatRaster, stars or sf object defining the prediction area (see
  Details)

- dist_space:

  "geographical", "feature" or "time". Should the distance be computed
  in geographic space, in the normalized multivariate predictor space or
  in temporal space? (see Details)

- CVtest:

  optional. list or vector. \#' @param cvfolds optional. list or vector.
  Either a list with the length of the number of cross-validation folds
  where each element contains the row indices of the data points used
  for testing during the cross validation iteration (i.e. held back
  data). Or a vector that contains the ID of the fold for each training
  point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm

- CVtrain:

  optional. A list, where each element contains the data points used for
  training during the cross validation iteration. Only required if
  CVtrain is not the opposite of CVtest. Relevant if some data points
  are excluded, e.g. when using
  [`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md).

- testdata:

  optional. object of class sf: Point data used for independent
  validation. May already include the predictor values if
  \`dist_space\`=feature.

- preddata:

  optional. object of class sf: Point data indicating the locations
  within the modeldomain to be used as target prediction points. Useful
  when the prediction objective is a subset of locations within the
  modeldomain rather than the whole area. May already include the
  predictor values if \`dist_space\`="feature".

- samplesize:

  numeric. How many prediction samples should be used?

- sampling:

  character. How to draw prediction samples? See
  [st_sample](https://r-spatial.github.io/sf/reference/st_sample.html)
  for modeldomains that are sf objects and
  [spatSample](https://rspatial.github.io/terra/reference/sample.html)
  for raster objects. Use sampling = "Fibonacci" for global applications
  (raster objects will be transformed to polygons in this case).

- variables:

  character vector defining the predictor variables used if
  dist_space="feature". If not provided all variables included in
  modeldomain are used.

- time_var:

  optional. character. Column that indicates the date. Only used if
  dist_space="time".

- time_unit:

  optional. Character. Unit for temporal distances See ?difftime.Only
  used if dist_space="time".

- algorithm:

  see [`knnx.dist`](https://rdrr.io/pkg/FNN/man/knn.dist.html) and
  [`knnx.index`](https://rdrr.io/pkg/FNN/man/knn.index.html)

- dist_fun:

  character. Currently covers \`euclidean\` (default), \`gower\`,
  \`mahalanobis\`, \`great_circle\` and \`abs_time\`. \`gower\` and
  \`mahalanobis\` only work with \`dist_space\`="feature", while
  \`great_circle\` only works with \`dist_space\`="geographical".
  \`mahalanobis\` takes into account correlation between predictor
  values. While \`euclidean\` and \`mahalanobis\` only work with
  numerical variables, \`gower\` also works with mixed data including
  numerical and categorical variables. For \`dist_space\`="time",
  currently only the absolute difference (\`abs_time\`) is implemented.
  For the geographical space, \`great_circle\` covers lon/lat
  coordinates, whereas \`euclidean\` only works with projected
  coordinates.

- scale_vars:

  boolean. Should variables be scaled? Only for
  \`dist_space\`="feature". Calculating Gower distances already includes
  scaling, and manually rescale the data is redundant. For other
  distances (Mahalanobis, Euclidean), scaling the data is important.
  Thus, TRUE by default.

- cvtrain:

  deprecated. Use \`CVtrain\` instead.

- cvfolds:

  deprecated. Use \`CVtest\` instead.

- type:

  deprecated. Use \`dist_space\` instead.

- timevar:

  deprecated. Use \`time_var\` instead.

## Value

A data.frame containing the distances. Unit of returned geographic
distances is meters. attributes contain W statistic between prediction
area and either sample data, CV folds or test data. See details.

## Details

The modeldomain is a sf polygon or a raster that defines the prediction
area. The function takes a regular point sample (amount defined by
samplesize) from the spatial extent (if no \`preddata\` are supplied).
If \`dist_space\` = "feature", the argument modeldomain has to be a
raster and include predictors. The only exception is when the provided
training data and preddata already include the predictor values. If not
provided they are extracted from the modeldomain raster. If some
predictors are categorical (i.e., of class factor or character), gower
distances will be used. W statistic describes the match between the
distributions. See Linnenbrink et al (2024) for further details.

## Note

See Meyer and Pebesma (2022) for an application of this plotting
function

## See also

[`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md)
[`knndm`](https://hannameyer.github.io/CAST/reference/knndm.md)

## Author

Hanna Meyer, Edzer Pebesma, Marvin Ludwig, Jan Linnenbrink

## Examples

``` r
if (FALSE) { # \dontrun{
library(CAST)
library(sf)
library(terra)
library(caret)
library(rnaturalearth)
library(ggplot2)

data(splotdata)
studyArea <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf")

########### Distance between training data and new data:
dist <- geodist(splotdata, studyArea)
# With density functions
plot(dist)
# Or ECDFs (relevant for nndm and knnmd methods)
plot(dist, stat="ecdf")

########### Distance between training data, new data and test data (here Chile):
plot(splotdata[,"Country"])
dist <- geodist(splotdata[splotdata$Country != "Chile",], studyArea,
                testdata = splotdata[splotdata$Country == "Chile",])
plot(dist)

########### Distance between training data, new data and CV folds:
folds <- createFolds(1:nrow(splotdata), k=3, returnTrain=FALSE)
dist <- geodist(x=splotdata, modeldomain=studyArea, CVtest=folds)
# Using density functions
plot(dist)
# Using ECDFs (relevant for nndm and knnmd methods)
plot(dist, stat="ecdf")

########### Distances in the feature space:
predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))
dist <- geodist(x = splotdata,
                modeldomain = predictors,
                dist_space = "feature",
                variables = c("bio_1","bio_12", "elev"))
plot(dist)

dist <- geodist(x = splotdata[splotdata$Country != "Chile",],
                modeldomain = predictors,
                testdata = splotdata[splotdata$Country == "Chile",],
                dist_space = "feature",
                variables=c("bio_1","bio_12", "elev"))
plot(dist)

############Distances in temporal space
library(lubridate)
library(ggplot2)
data(cookfarm)
dat <- st_as_sf(cookfarm,coords=c("Easting","Northing"))
st_crs(dat) <- 26911
trainDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2010,]
predictionDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2011,]
trainDat$week <- lubridate::week(trainDat$Date)
CVtest <- CreateSpacetimeFolds(trainDat,time_var = "week")

dist <- geodist(trainDat,preddata = predictionDat,CVtest = CVtest$indexOut,
   dist_space="time",time_unit="days")
plot(dist)+ xlim(0,10)


############ Example for a random global dataset
############ (refer to figure in Meyer and Pebesma 2022)

### Define prediction area (here: global):
ee <- st_crs("+proj=eqearth")
co <- ne_countries(returnclass = "sf")
co.ee <- st_transform(co, ee)

### Simulate a spatial random sample
### (alternatively replace pts_random by a real sampling dataset (see Meyer and Pebesma 2022):
sf_use_s2(FALSE)
pts_random <- st_sample(co.ee, 2000, exact=FALSE)

### See points on the map:
ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_random, color = "#F8766D",size=0.5, shape=3) +
  guides(fill = "none", col = "none") +
  labs(x = NULL, y = NULL)

### plot distances:
dist <- geodist(pts_random,co.ee)
plot(dist) + scale_x_log10(labels=round)




} # }
```
