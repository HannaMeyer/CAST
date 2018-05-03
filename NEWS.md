# `CAST` 0.2.0

* bugfix: `ffs` with option withinSE=TRUE did not choose a model as "best model" if it was within the SE of a model that was trained in an earlier run but had the same number of variables. This bug is fixed and if withinSE=TRUE ffs now only compares the performance to models that use less variables (e.g. if a model using 5 variables is better than a model using 4 variables but still in the SE of the 4-variable model, then the 4-variable model is rated as the better model).

* new feature: `plot_ffs` plots the results of ffs to visualize how the performance changes according to model run and the number of variables being used.

# `CAST` 0.1.0

Initial public version on CRAN

