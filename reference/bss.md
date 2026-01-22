# Best subset feature selection

Evaluate all combinations of predictors during model training

## Usage

``` r
bss(
  predictors,
  response,
  method = "rf",
  metric = ifelse(is.factor(response), "Accuracy", "RMSE"),
  maximize = ifelse(metric == "RMSE", FALSE, TRUE),
  globalval = FALSE,
  trControl = caret::trainControl(),
  tuneLength = 3,
  tuneGrid = NULL,
  seed = 100,
  verbose = TRUE,
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

- trControl:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- tuneLength:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- tuneGrid:

  see [`train`](https://rdrr.io/pkg/caret/man/train.html)

- seed:

  A random number

- verbose:

  Logical. Should information about the progress be printed?

- ...:

  arguments passed to the classification or regression routine (such as
  randomForest).

## Value

A list of class train. Beside of the usual train content the object
contains the vector "selectedvars" and "selectedvars_perf" that give the
best variables selected as well as their corresponding performance. It
also contains "perf_all" that gives the performance of all model runs.

## Details

bss is an alternative to
[`ffs`](https://hannameyer.github.io/CAST/reference/ffs.md) and ideal if
the training set is small. Models are iteratively fitted using all
different combinations of predictor variables. Hence, 2^X models are
calculated. Don't try running bss on very large datasets because the
computation time is much higher compared to
[`ffs`](https://hannameyer.github.io/CAST/reference/ffs.md).

The internal cross validation can be run in parallel. See information on
parallel processing of carets train functions for details.

## Note

This variable selection is particularly suitable for spatial cross
validations where variable selection MUST be based on the performance of
the model for predicting new spatial units. Note that bss is very slow
since all combinations of variables are tested. A more time efficient
alternative is the forward feature selection
([`ffs`](https://hannameyer.github.io/CAST/reference/ffs.md)).

## See also

[`train`](https://rdrr.io/pkg/caret/man/train.html),[`ffs`](https://hannameyer.github.io/CAST/reference/ffs.md),
[`trainControl`](https://rdrr.io/pkg/caret/man/trainControl.html),[`CreateSpacetimeFolds`](https://hannameyer.github.io/CAST/reference/CreateSpacetimeFolds.md),
[`nndm`](https://hannameyer.github.io/CAST/reference/nndm.md)

## Author

Hanna Meyer

## Examples

``` r
if (FALSE) { # \dontrun{
data(iris)
bssmodel <- bss(iris[,1:4],iris$Species)
bssmodel$perf_all
plot(bssmodel)
} # }
```
