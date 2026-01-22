# Calculate Dissimilarity Index of training data

This function estimates the Dissimilarity Index (DI) within the training
data set used for a prediction model. Optionally, the local point
density can also be calculated. Predictors can be weighted based on the
internal variable importance of the machine learning algorithm used for
model training.

## Usage

``` r
trainDI(
  model = NA,
  train = NULL,
  variables = "all",
  weight = NA,
  CVtest = NULL,
  CVtrain = NULL,
  method = "L2",
  useWeight = TRUE,
  useCV = TRUE,
  LPD = FALSE,
  verbose = TRUE,
  algorithm = "brute"
)
```

## Arguments

- model:

  A train object created with caret used to extract weights from (based
  on variable importance) as well as cross-validation folds

- train:

  A data.frame containing the data used for model training. Only
  required when no model is given

- variables:

  character vector of predictor variables. if "all" then all variables
  of the model are used or if no model is given then of the train
  dataset.

- weight:

  A data.frame containing weights for each variable. Only required if no
  model is given.

- CVtest:

  list or vector. Either a list where each element contains the data
  points used for testing during the cross validation iteration (i.e.
  held back data). Or a vector that contains the ID of the fold for each
  training point. Only required if no model is given.

- CVtrain:

  list. Each element contains the data points used for training during
  the cross validation iteration (i.e. held back data). Only required if
  no model is given and only required if CVtrain is not the opposite of
  CVtest (i.e. if a data point is not used for testing, it is used for
  training). Relevant if some data points are excluded, e.g. when using
  [`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md).

- method:

  Character. Method used for distance calculation. Currently euclidean
  distance (L2) and Mahalanobis distance (MD) are implemented but only
  L2 is tested. Note that MD takes considerably longer.

- useWeight:

  Logical. Only if a model is given. Weight variables according to
  importance in the model?

- useCV:

  Logical. Only if a model is given. Use the CV folds to calculate the
  DI threshold?

- LPD:

  Logical. Indicates whether the local point density should be
  calculated or not.

- verbose:

  Logical. Print progress or not?

- algorithm:

  see [`knnx.dist`](https://rdrr.io/pkg/FNN/man/knn.dist.html) and
  [`knnx.index`](https://rdrr.io/pkg/FNN/man/knn.index.html)

## Value

A list of class `trainDI` containing:

- train:

  A data frame containing the training data

- weight:

  A data frame with weights based on the variable importance.

- variables:

  Names of the used variables

- catvars:

  Which variables are categorial

- scaleparam:

  Scaling parameters. Output from `scale`

- trainDist_avrg:

  A data frame with the average distance of each training point to every
  other point

- trainDist_avrgmean:

  The mean of trainDist_avrg. Used for normalizing the DI

- trainDI:

  Dissimilarity Index of the training data

- threshold:

  The DI threshold used for inside/outside AOA

- trainLPD:

  LPD of the training data

- avrgLPD:

  Average LPD of the training data

## Note

This function is called within
[`aoa`](https://hannameyer.github.io/CAST/reference/aoa.md) to estimate
the DI and AOA of new data. However, it may also be used on its own if
only the DI of training data is of interest, or to facilitate a
parallelization of
[`aoa`](https://hannameyer.github.io/CAST/reference/aoa.md) by avoiding
a repeated calculation of the DI within the training data.

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
library(sf)
library(terra)
library(caret)
library(CAST)

# prepare sample data:
data("splotdata")
splotdata = st_drop_geometry(splotdata)

# train a model:
set.seed(100)
model <- caret::train(splotdata[,6:16],
                      splotdata$Species_richness,
                      importance=TRUE, tuneLength=1, ntree = 15, method = "rf",
                      trControl = trainControl(method="cv", number=5, savePredictions=T))
# variable importance is used for scaling predictors
plot(varImp(model,scale=FALSE))

# calculate the DI of the trained model:
DI = trainDI(model=model)
plot(DI)

#...or calculate the DI and LPD of the trained model:
# DI = trainDI(model=model, LPD = TRUE)

# the DI can now be used to compute the AOA (here with LPD):
studyArea = rast(system.file("extdata/predictors_chile.tif", package = "CAST"))
AOA = aoa(studyArea, model = model, trainDI = DI, LPD = TRUE, maxLPD = 1)
print(AOA)
plot(AOA)
plot(AOA$AOA)
plot(AOA$LPD)
} # }
```
