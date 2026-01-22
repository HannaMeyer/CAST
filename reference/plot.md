# Plot CAST classes

Generic plot function for CAST Classes

A plotting function for a forward feature selection result. Each point
is the mean performance of a model run. Error bars represent the
standard errors from cross validation. Marked points show the best model
from each number of variables until a further variable could not improve
the results. If type=="selected", the contribution of the selected
variables to the model performance is shown.

Density plot of nearest neighbor distances in geographic space or
feature space between training data as well as between training data and
prediction locations. Optional, the nearest neighbor distances between
training data and test data or between training data and CV iterations
is shown. The plot can be used to check the suitability of a chosen CV
method to be representative to estimate map accuracy.

Plot the DI/LPD and errormetric from Cross-Validation with the modeled
relationship

## Usage

``` r
# S3 method for class 'trainDI'
plot(x, ...)

# S3 method for class 'aoa'
plot(x, samplesize = 1000, variable = "DI", ...)

# S3 method for class 'nndm'
plot(x, type = "strict", stat = "ecdf", ...)

# S3 method for class 'knndm'
plot(x, type = "strict", stat = "ecdf", ...)

# S3 method for class 'ffs'
plot(
  x,
  plotType = "all",
  palette = rainbow,
  reverse = FALSE,
  marker = "black",
  size = 1.5,
  lwd = 0.5,
  pch = 21,
  ...
)

# S3 method for class 'geodist'
plot(x, unit = "m", stat = "density", ...)

# S3 method for class 'errorModel'
plot(x, ...)
```

## Arguments

- x:

  errorModel, see
  [`DItoErrormetric`](https://hannameyer.github.io/CAST/reference/errorProfiles.md)

- ...:

  other params

- samplesize:

  numeric. How many prediction samples should be plotted?

- variable:

  character. Variable for which to generate the density plot. 'DI' or
  'LPD'

- type:

  String, defaults to "strict" to show the original nearest neighbour
  distance definitions in the legend. Alternatively, set to "simple" to
  have more intuitive labels.

- stat:

  "density" for density plot or "ecdf" for empirical cumulative
  distribution function plot.

- plotType:

  character. Either "all" or "selected"

- palette:

  A color palette

- reverse:

  Character. Should the palette be reversed?

- marker:

  Character. Color to mark the best models

- size:

  Numeric. Size of the points

- lwd:

  Numeric. Width of the error bars

- pch:

  Numeric. Type of point marking the best models

- unit:

  character. Only if type=="geo" and only applied to the plot.
  Supported: "m" or "km".

## Value

a ggplot

a ggplot

## Author

Marvin Ludwig, Hanna Meyer

Carles Milà

## Examples

``` r
if (FALSE) { # \dontrun{
data(splotdata)
splotdata <- st_drop_geometry(splotdata)
ffsmodel <- ffs(splotdata[,6:16], splotdata$Species_richness, ntree = 10)
plot(ffsmodel)
#plot performance of selected variables only:
plot(ffsmodel,plotType="selected")
} # }
```
