# Forward feature selection

A simple forward feature selection algorithm

## Usage

``` r
ffs(
  predictors,
  response,
  method = "rf",
  metric = ifelse(is.factor(response), "Accuracy", "RMSE"),
  maximize = ifelse(metric == "RMSE", FALSE, TRUE),
  globalval = FALSE,
  withinSE = FALSE,
  minVar = 2,
  trControl = caret::trainControl(),
  tuneLength = 3,
  tuneGrid = NULL,
  seed = sample(1:1000, 1),
  verbose = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- predictors:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- response:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- method:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- metric:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- maximize:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- globalval:

  Logical. Should models be evaluated based on 'global' performance? See
  [`global_validation`](https://hannameyer.github.io/CAST/reference/global_validation.md)

- withinSE:

  Logical Models are only selected if they are better than the currently
  best models Standard error

- minVar:

  Numeric. Number of variables to combine for the first selection. See
  Details.

- trControl:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- tuneLength:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- tuneGrid:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- seed:

  A random number used for model training

- verbose:

  Logical. Should information about the progress be printed?

- cores:

  Numeric. If \> 2, mclapply will be used. see
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html)

- ...:

  arguments passed to the classification or regression routine (such as
  randomForest).

## Value

A list of class train. Beside of the usual train content the object
contains the vector "selectedvars" and "selectedvars_perf" that give the
order of the best variables selected as well as their corresponding
performance (starting from the first two variables). It also contains
"perf_all" that gives the performance of all model runs.

## Details

Models with two predictors are first trained using all possible pairs of
predictor variables. The best model of these initial models is kept. On
the basis of this best model the predictor variables are iteratively
increased and each of the remaining variables is tested for its
improvement of the currently best model. The process stops if none of
the remaining variables increases the model performance when added to
the current best model.

The forward feature selection can be run in parallel with forking on
Linux systems (mclapply). Each fork computes a model, which drastically
speeds up the runtime - especially of the initial predictor search. The
internal cross validation can be run in parallel on all systems. See
information on parallel processing of carets train functions for
details.

Using withinSE will favour models with less variables and probably
shorten the calculation time

Per Default, the ffs starts with all possible 2-pair combinations.
minVar allows to start the selection with more than 2 variables, e.g.
minVar=3 starts the ffs testing all combinations of 3 (instead of 2)
variables first and then increasing the number. This is important for
e.g. neural networks that often cannot make sense of only two variables.
It is also relevant if it is assumed that the optimal variables can only
be found if more than 2 are considered at the same time.

## Note

This variable selection is particularly suitable for spatial cross
validations where variable selection MUST be based on the performance of
the model for predicting new spatial units. See Meyer et al. (2018) and
Meyer et al. (2019) for further details.

## References

- Gasch, C.K., Hengl, T., Gräler, B., Meyer, H., Magney, T., Brown, D.J.
  (2015): Spatio-temporal interpolation of soil water, temperature, and
  electrical conductivity in 3D+T: the Cook Agronomy Farm data set.
  Spatial Statistics 14: 70-90.

- Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018):
  Improving performance of spatio-temporal machine learning models using
  forward feature selection and target-oriented validation.
  Environmental Modelling & Software 101: 1-9.
  [doi:10.1016/j.envsoft.2017.12.001](https://doi.org/10.1016/j.envsoft.2017.12.001)

- Meyer, H., Reudenbach, C., Wöllauer, S., Nauss, T. (2019): Importance
  of spatial predictor variable selection in machine learning
  applications - Moving from data reproduction to spatial prediction.
  Ecological Modelling. 411, 108815.
  [doi:10.1016/j.ecolmodel.2019.108815](https://doi.org/10.1016/j.ecolmodel.2019.108815)
  .

- Ludwig, M., Moreno-Martinez, A., Hölzel, N., Pebesma, E., Meyer, H.
  (2023): Assessing and improving the transferability of current global
  spatial prediction models. Global Ecology and Biogeography.
  [doi:10.1111/geb.13635](https://doi.org/10.1111/geb.13635) .

## See also

[`train`](https://rdrr.io/pkg/caret/man/train.html),[`bss`](https://hannameyer.github.io/CAST/reference/bss.md),
[`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html),[`CreateSpacetimeFolds`](https://hannameyer.github.io/CAST/reference/CreateSpacetimeFolds.md),[`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md)

## Author

Hanna Meyer

## Examples

``` r
if (FALSE) { # \dontrun{
data(splotdata)
ffsmodel <- ffs(splotdata[,6:12], splotdata$Species_richness, ntree = 20)

ffsmodel$selectedvars
ffsmodel$selectedvars_perf
plot(ffsmodel)
#or only selected variables:
plot(ffsmodel,plotType="selected")
} # }

# or perform model with target-oriented validation (LLO CV)
#the example is described in Gasch et al. (2015). The ffs approach for this dataset is described in
#Meyer et al. (2018). Due to high computation time needed, only a small and thus not robust example
#is shown here.

if (FALSE) { # \dontrun{
# run the model on three cores (see vignette for details):
library(doParallel)
library(lubridate)
cl <- makeCluster(3)
registerDoParallel(cl)

#load and prepare dataset:
data(cookfarm)
trainDat <- cookfarm[cookfarm$altitude==-0.3&
  year(cookfarm$Date)==2012&week(cookfarm$Date)%in%c(13:14),]

#visualize dataset:
ggplot(data = trainDat, aes(x=Date, y=VW)) + geom_line(aes(colour=SOURCEID))

#create folds for Leave Location Out Cross Validation:
set.seed(10)
indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID",k=3)
ctrl <- trainControl(method="cv",index = indices$index)

#define potential predictors:
predictors <- c("DEM","TWI","BLD","Precip_cum","cday","MaxT_wrcc",
"Precip_wrcc","NDRE.M","Bt","MinT_wrcc","Northing","Easting")

#run ffs model with Leave Location out CV
set.seed(10)
ffsmodel <- ffs(trainDat[,predictors],trainDat$VW,method="rf",
tuneLength=1,trControl=ctrl)
ffsmodel
plot(ffsmodel)
#or only selected variables:
plot(ffsmodel,plotType="selected")

#compare to model without ffs:
model <- train(trainDat[,predictors],trainDat$VW,method="rf",
tuneLength=1, trControl=ctrl)
model
stopCluster(cl)
} # }

if (FALSE) { # \dontrun{
## on linux machines, you can also run the ffs in parallel with forks:
data("splotdata")
spatial_cv = CreateSpacetimeFolds(splotdata, spacevar = "Biome", k = 5)
ctrl <- trainControl(method="cv",index = spatial_cv$index)

ffsmodel <- ffs(predictors = splotdata[,6:16],
               response = splotdata$Species_richness,
               tuneLength = 1,
               method = "rf",
               trControl = ctrl,
               ntree = 20,
               seed = 1,
               cores = 4)
} # }

```
