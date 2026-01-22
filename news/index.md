# Changelog

## `CAST` 1.0.4

- bug fix: formatting issues in tests
- references updated

## `CAST` 1.0.3

CRAN release: 2025-01-09

- bug fix: default algorithm for FNN functions changed

## `CAST` 1.0.2

CRAN release: 2024-06-14

- bug fix: tests run conditionally

## `CAST` 1.0.1

CRAN release: 2024-04-25

- bug fix: fix failed tests in global_validation

## `CAST` 1.0.0

CRAN release: 2024-04-08

- new features:
  - calculate local point density within AOA
  - option of spatial error profiles (errorProfiles with
    variable=“geodist”)
  - normalize_DI for a more intuitive interpretation
  - geodist allows calculating temporal distances
  - ffs now can be run in parallel (Linux only)
  - vignette “Cross-validation methods in CAST”
  - knndm in feature space (experimental)
  - nndm in feature space (experimental)
- modifications:
  - function DItoErrormetric renamed to errorProfiles and allows for
    other dissimilarity measures
  - Improvement and homogenization of plotting methods for nndm, knndm
    and geodist objects
  - aoa and trainDI `weight` now allows list input
  - vignette on Introduction to CAST updated
- deprecated: *plot_geodist (replaced by plot.geodist)* plot_ffs
  (replaced by plot.ffs) \*calibrate_aoa (replaced by errorProfiles)

## `CAST` 0.9.0

CRAN release: 2024-01-09

- new features:
  - CAST functions now return classes with generic plotting and printing
  - new dataset for examples, tutorials and testing: data(splotdata)
- modifications:
  - calibrate_aoa is now DItoErrormetric and returns a model (see
    function documentation)
  - plot_geodist is now geodist. The result can be visualized with
    plot()
  - plot_ffs is now plot(ffs)
- bug fix:
  - fix issue [\#65](https://github.com/HannaMeyer/CAST/issues/65)
    (threshold)
- deprecated (soon):
  - plot_geodist, plot_ffs, calibrate_aoa

## `CAST` 0.8.1

CRAN release: 2023-05-30

- bugfix:
  - failed checks on Fedora 34 fixed

## `CAST` 0.8.0

CRAN release: 2023-05-21

- new features:
  - knndm as an alternative to nndm for large training data
- modifications:
  - transition from raster to terra

## `CAST` 0.7.1

CRAN release: 2023-02-04

- new features:
  - Mahalanobis distance for AOA assessment as option
- modifications:
  - faster estimation of the AOA
  - parallel option for AOA deprecated (see vignette)
    - delineation of the default threshold fixed as suggested in
      github.com/HannaMeyer/CAST/issues/46
- bugfix:
  - fixed issue github.com/ropensci/rnaturalearth/issues/69

## `CAST` 0.7.0

CRAN release: 2022-08-24

- new feature:
  - nndm cross-validation as suggested by Milà et al. (2022)
- modifications
  - plot_geodist works with NNDM
  - trainDI works with NNDM
  - rename of parameter folds in AOA and trainDI

## `CAST` 0.6.0

CRAN release: 2022-03-17

- new feature:
  - trainDI allows to calculate the DI of the training dataset
    separately from the aoa function
  - plot and print functions for the AOA
  - function to plot nearest neighbor distance distributions in
    geographic and feature space
  - function global_validation added
- modifications
  - extensive restructuring of the AOA function
  - ffs and bss can be used with global_validation
- bugfix:
  - error in manual assignment of weights fixed

## `CAST` 0.5.1

CRAN release: 2021-04-07

- resolved dependence on package “GSIF” which was removed from the CRAN
  repository

## `CAST` 0.5.0

CRAN release: 2021-02-19

- new feature:
  - AOA can run in parallel
  - calibration of the DI (calibrate_aoa)
- bugfix:
  - aoa will work now with large training sets
- modifications:
  - default threshold of AOA changed

## `CAST` 0.4.2

CRAN release: 2020-07-17

- new feature:
  - aoa now working with categorical variables
- bugfix:
  - fixed error in ffs when \>170 variables are used
- minor changes:
  - changed order of parameters in aoa
  - tutorial “Introduction to CAST” improved

## `CAST` 0.4.1

CRAN release: 2020-05-19

- new feature:
  - vignette: tutorial introducing the “area of applicability”
  - variable threshold for aoa
  - various modifications in aoa in line with submitted paper

## `CAST` 0.4.0

CRAN release: 2020-04-06

- new feature:
  - new function “aoa”: quantify and visualize the area of applicability
    of spatial prediction models
  - “minVar” in ffs: Instead of always starting with 2-pair
    combinations, ffs can now also be started with combinations of more
    variables (e.g starting with all combinations of 3)
- bugfix:
  - ffs failed for “svmLinear” in previous version because of S4 class
    issues. Fixed now.

## `CAST` 0.3.1

CRAN release: 2018-11-19

- bugfix:
  - CreateSpaceTimeFolds accepts tibbles
  - CreateSpaceTimeFolds automatically reduces k if necessary
  - ffs accepts further arguments taken by caret::train
- new feature: plot_ffs has option to plot selected variables only

## `CAST` 0.3.0

CRAN release: 2018-10-11

- new feature: Best subset selection (bss) with target-oriented
  validation as (very slow but very reliable) alternative to ffs

- minor adaptations: verbose option included, improved examples for ffs

- bugfix: minor adaptations done for usage with plsr

## `CAST` 0.2.1

CRAN release: 2018-07-12

- new feature: Introduction to CAST is included as a vignette.

- bugfix: minor error fixed in using user defined metrics for model
  selection.

## `CAST` 0.2.0

CRAN release: 2018-05-03

- bugfix: `ffs` with option withinSE=TRUE did not choose a model as
  “best model” if it was within the SE of a model that was trained in an
  earlier run but had the same number of variables. This bug is fixed
  and if withinSE=TRUE ffs now only compares the performance to models
  that use less variables (e.g. if a model using 5 variables is better
  than a model using 4 variables but still in the SE of the 4-variable
  model, then the 4-variable model is rated as the better model).

- new feature: `plot_ffs` plots the results of ffs to visualize how the
  performance changes according to model run and the number of variables
  being used.

## `CAST` 0.1.0

CRAN release: 2018-01-09

Initial public version on CRAN
