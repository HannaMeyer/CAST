# Calculate euclidean nearest neighbor distances in geographic space or feature space

Calculates nearest neighbor distances in geographic space or feature
space between training data as well as between training data and
prediction locations. Optional, the nearest neighbor distances between
training data and test data or between training data and CV iterations
is computed.

## Usage

``` r
geodist(
  x,
  modeldomain = NULL,
  type = "geo",
  cvfolds = NULL,
  cvtrain = NULL,
  testdata = NULL,
  preddata = NULL,
  samplesize = 2000,
  sampling = "regular",
  variables = NULL,
  timevar = NULL,
  time_unit = "auto",
  algorithm = "brute"
)
```

## Arguments

- x:

  object of class sf, training data locations

- modeldomain:

  SpatRaster, stars or sf object defining the prediction area (see
  Details)

- type:

  "geo" or "feature". Should the distance be computed in geographic
  space or in the normalized multivariate predictor space (see Details)

- cvfolds:

  optional. list or vector. Either a list where each element contains
  the data points used for testing during the cross validation iteration
  (i.e. held back data). Or a vector that contains the ID of the fold
  for each training point. See e.g. ?createFolds or
  ?CreateSpacetimeFolds or ?nndm

- cvtrain:

  optional. List of row indices of x to fit the model to in each CV
  iteration. If cvtrain is null but cvfolds is not, all samples but
  those included in cvfolds are used as training data

- testdata:

  optional. object of class sf: Point data used for independent
  validation

- preddata:

  optional. object of class sf: Point data indicating the locations
  within the modeldomain to be used as target prediction points. Useful
  when the prediction objective is a subset of locations within the
  modeldomain rather than the whole area.

- samplesize:

  numeric. How many prediction samples should be used?

- sampling:

  character. How to draw prediction samples? See
  [spsample](https://edzer.github.io/sp/reference/spsample.html). Use
  sampling = "Fibonacci" for global applications.

- variables:

  character vector defining the predictor variables used if
  type="feature. If not provided all variables included in modeldomain
  are used.

- timevar:

  optional. character. Column that indicates the date. Only used if
  type="time".

- time_unit:

  optional. Character. Unit for temporal distances See ?difftime.Only
  used if type="time".

- algorithm:

  see [`knnx.dist`](https://rdrr.io/pkg/FNN/man/knn.dist.html) and
  [`knnx.index`](https://rdrr.io/pkg/FNN/man/knn.index.html)

## Value

A data.frame containing the distances. Unit of returned geographic
distances is meters. attributes contain W statistic between prediction
area and either sample data, CV folds or test data. See details.

## Details

The modeldomain is a sf polygon or a raster that defines the prediction
area. The function takes a regular point sample (amount defined by
samplesize) from the spatial extent. If type = "feature", the argument
modeldomain (and if provided then also the testdata and/or preddata) has
to include predictors. Predictor values for x, testdata and preddata are
optional if modeldomain is a raster. If not provided they are extracted
from the modeldomain rasterStack. If some predictors are categorical
(i.e., of class factor or character), gower distances will be used. W
statistic describes the match between the distributions. See Linnenbrink
et al (2023) for further details.

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
dist <- geodist(x=splotdata, modeldomain=studyArea, cvfolds=folds)
# Using density functions
plot(dist)
# Using ECDFs (relevant for nndm and knnmd methods)
plot(dist, stat="ecdf")

########### Distances in the feature space:
predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))
dist <- geodist(x = splotdata,
                modeldomain = predictors,
                type = "feature",
                variables = c("bio_1","bio_12", "elev"))
plot(dist)

dist <- geodist(x = splotdata[splotdata$Country != "Chile",],
                modeldomain = predictors, cvfolds = folds,
                testdata = splotdata[splotdata$Country == "Chile",],
                type = "feature",
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
cvfolds <- CreateSpacetimeFolds(trainDat,timevar = "week")

dist <- geodist(trainDat,preddata = predictionDat,cvfolds = cvfolds$indexOut,
   type="time",time_unit="days")
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
