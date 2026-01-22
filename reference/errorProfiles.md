# Model and inspect the relationship between the prediction error and measures of dissimilarities and distances

Performance metrics are calculated for moving windows of dissimilarity
values based on cross-validated training data

## Usage

``` r
errorProfiles(
  model,
  trainDI = NULL,
  locations = NULL,
  variable = "DI",
  multiCV = FALSE,
  length.out = 10,
  window.size = 5,
  calib = "scam",
  method = "L2",
  useWeight = TRUE,
  k = 6,
  m = 2
)
```

## Arguments

- model:

  the model used to get the AOA

- trainDI:

  the result of
  [`trainDI`](https://hannameyer.github.io/CAST/reference/trainDI.md) or
  aoa object [`aoa`](https://hannameyer.github.io/CAST/reference/aoa.md)

- locations:

  Optional. sf object for the training data used in model. Only used if
  variable=="geodist". Note that they must be in the same order as
  model\$trainingData.

- variable:

  Character. Which dissimilarity or distance measure to use for the
  error metric. Current options are "DI" or "LPD"

- multiCV:

  Logical. Re-run model fitting and validation with different CV
  strategies. See details.

- length.out:

  Numeric. Only used if multiCV=TRUE. Number of cross-validation folds.
  See details.

- window.size:

  Numeric. Size of the moving window. See
  [`rollapply`](https://rdrr.io/pkg/zoo/man/rollapply.html).

- calib:

  Character. Function to model the DI/LPD~performance relationship.
  Currently lm and scam are supported

- method:

  Character. Method used for distance calculation. Currently euclidean
  distance (L2) and Mahalanobis distance (MD) are implemented but only
  L2 is tested. Note that MD takes considerably longer. See ?aoa for
  further explanation

- useWeight:

  Logical. Only if a model is given. Weight variables according to
  importance in the model?

- k:

  Numeric. See mgcv::s

- m:

  Numeric. See mgcv::s

## Value

A scam, linear model or exponential model

## Details

If multiCV=TRUE the model is re-fitted and validated by length.out new
cross-validations where the cross-validation folds are defined by
clusters in the predictor space, ranging from three clusters to LOOCV.
Hence, a large range of dissimilarity values is created during
cross-validation. If the AOA threshold based on the calibration data
from multiple CV is larger than the original AOA threshold (which is
likely if extrapolation situations are created during CV), the AOA
threshold changes accordingly. See Meyer and Pebesma (2021) for the full
documentation of the methodology.

## References

Meyer, H., Pebesma, E. (2021): Predicting into unknown space? Estimating
the area of applicability of spatial prediction models.
[doi:10.1111/2041-210X.13650](https://doi.org/10.1111/2041-210X.13650)

## See also

[`aoa`](https://hannameyer.github.io/CAST/reference/aoa.md)

## Author

Hanna Meyer, Marvin Ludwig, Fabian Schumacher

## Examples

``` r
if (FALSE) { # \dontrun{
library(CAST)
library(sf)
library(terra)
library(caret)

data(splotdata)
predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

model <- caret::train(st_drop_geometry(splotdata)[,6:16], splotdata$Species_richness,
   ntree = 10, trControl = trainControl(method = "cv", savePredictions = TRUE))

AOA <- aoa(predictors, model, LPD = TRUE, maxLPD = 1)

### DI ~ error
errormodel_DI <- errorProfiles(model, AOA, variable = "DI")
plot(errormodel_DI)
summary(errormodel_DI)

expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
plot(expected_error_DI)

### LPD ~ error
errormodel_LPD <- errorProfiles(model, AOA, variable = "LPD")
plot(errormodel_LPD)
summary(errormodel_DI)

expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)
plot(expected_error_LPD)

### geodist ~ error
errormodel_geodist = errorProfiles(model, locations=splotdata, variable = "geodist")
plot(errormodel_geodist)
summary(errormodel_DI)

dist <- terra::distance(predictors[[1]],vect(splotdata))
names(dist) <- "geodist"
expected_error_DI <- terra::predict(dist, errormodel_geodist)
plot(expected_error_DI)


### with multiCV = TRUE (for DI ~ error)
errormodel_DI = errorProfiles(model, AOA, multiCV = TRUE, length.out = 3, variable = "DI")
plot(errormodel_DI)

expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
plot(expected_error_DI)

# mask AOA based on new threshold from multiCV
mask_aoa = terra::mask(expected_error_DI, AOA$DI > attr(errormodel_DI, 'AOA_threshold'),
  maskvalues = 1)
plot(mask_aoa)
} # }

```
