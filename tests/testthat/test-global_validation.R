test_that("global_validation correctly handles missing predictions", {

  data("iris")
  set.seed(123)
  ctrl <- caret::trainControl(method="cv")
  model <- caret::train(iris[,c("Sepal.Width", "Petal.Length", "Petal.Width")],
                        iris[,c("Sepal.Length")],
                        method="rf", trControl=ctrl, ntree=10)
  expect_error(global_validation(model))
})

test_that("global_validation works with caret regression", {

  data("iris")
  set.seed(123)
  ctrl <- caret::trainControl(method="cv", savePredictions="final")
  model <- caret::train(iris[,c("Sepal.Width", "Petal.Length", "Petal.Width")],
                        iris[,c("Sepal.Length")],
                        method="rf", trControl=ctrl, ntree=10)
  expect_equal(global_validation(model),
               c("RMSE"=0.3307870, "Rsquared"=0.8400544, "MAE"=0.2621827))

})

test_that("global_validation works with caret classification", {

  data("iris")
  set.seed(123)
  ctrl <- caret::trainControl(method="cv", savePredictions="final")
  model <- caret::train(iris[,c("Sepal.Width", "Petal.Length", "Petal.Width", "Sepal.Length")],
                        iris[,c("Species")],
                        method="rf", trControl=ctrl, ntree=10)
  expect_equal(global_validation(model)[1:2],
               c("Accuracy"=0.96, "Kappa"=0.94))

})

test_that("global_validation works with CreateSpacetimeFolds", {

  data("iris")
  set.seed(123)
  iris$folds <- sample(rep(1:10, ceiling(nrow(iris)/10)), nrow(iris))
  indices <- CreateSpacetimeFolds(iris, "folds")
  ctrl <- caret::trainControl(method="cv", savePredictions="final", index = indices$index)
  model <- caret::train(iris[,c("Sepal.Width", "Petal.Length", "Petal.Width", "Sepal.Length")],
                        iris[,c("Species")],
                        method="rf", trControl=ctrl, ntree=10)
  expect_equal(global_validation(model)[1:2],
               c("Accuracy"=0.96, "Kappa"=0.94))
})
