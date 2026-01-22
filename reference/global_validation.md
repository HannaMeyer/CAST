# Evaluate 'global' cross-validation

Calculate validation metric using all held back predictions at once

## Usage

``` r
global_validation(model)
```

## Arguments

- model:

  an object of class [`train`](https://rdrr.io/pkg/caret/man/train.html)

## Value

regression
([`postResample`](https://rdrr.io/pkg/caret/man/postResample.html)) or
classification
([`confusionMatrix`](https://rdrr.io/pkg/caret/man/confusionMatrix.html))
statistics

## Details

Relevant when folds are not representative for the entire area of
interest. In this case, metrics like R2 are not meaningful since it
doesn't reflect the general ability of the model to explain the entire
gradient of the response. Comparable to LOOCV, predictions from all held
back folds are used here together to calculate validation statistics.

## See also

[`CreateSpacetimeFolds`](https://hannameyer.github.io/CAST/reference/CreateSpacetimeFolds.md)

## Author

Hanna Meyer

## Examples

``` r
if (FALSE) { # \dontrun{
library(caret)
data(cookfarm)
dat <- cookfarm[sample(1:nrow(cookfarm),500),]
indices <- CreateSpacetimeFolds(dat,"SOURCEID","Date")
ctrl <- caret::trainControl(method="cv",index = indices$index,savePredictions="final")
model <- caret::train(dat[,c("DEM","TWI","BLD")],dat$VW, method="rf", trControl=ctrl, ntree=10)
global_validation(model)
} # }
```
