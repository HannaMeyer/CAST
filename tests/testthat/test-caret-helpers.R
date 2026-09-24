skip_if_not_installed("caret")
source("data-fixture.R")

test_that("caret helpers work with caret models from fixture", {
  data <- loaddata()
  model <- data$model
  variables <- data$variables
  # .caret_get_data
  df <- .caret_get_data(model)
  expect_s3_class(df, "data.frame")
  expect_true(nrow(df) > 0)
  # .caret_get_variables
  vars <- .caret_get_variables(model)
  expect_type(vars, "character")
  expect_false(any(vars == ".outcome"))
  expect_true(length(vars) > 0)
  # .caret_get_folds
  folds <- .caret_get_folds(model)
  expect_type(folds, "list")
  expect_true(all(c("CVtrain", "CVtest") %in% names(folds)))
  if (!is.null(model$control) && !is.null(model$control$method) && tolower(model$control$method) == "cv") {
    expect_true(!is.null(folds$CVtrain) || !is.null(folds$CVtest))
  } else {
    expect_null(folds$CVtrain)
    expect_null(folds$CVtest)
  }
  # .caret_get_weights
  weights <- .caret_get_weights(model, vars)
  expect_s3_class(weights, "data.frame")
  expect_equal(ncol(weights), length(vars))
  expect_true(all(names(weights) %in% vars))
})
