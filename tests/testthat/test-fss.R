
test_that("ffs works with default arguments and the splotopen dataset (numerical only)",{
  data("splotdata")
  splotdata = splotdata |> sf::st_drop_geometry()
  set.seed(1)
  selection = ffs(predictors = splotdata[,6:16],
                  response = splotdata$Species_richness,
                  seed = 1,
                  verbose = FALSE,
                  ntree = 5)


  expect_identical(selection$selectedvars, c("bio_6", "bio_12", "bio_15"))
  expect_identical(selection$metric, "RMSE")
  expect_identical(selection$maximize, FALSE)

})




test_that("ffs works with default arguments and the splotopen dataset (include categorial)",{
  data("splotdata")
  splotdata = splotdata |> sf::st_drop_geometry()
  set.seed(1)
  selection = ffs(predictors = splotdata[,c(4,6:16)],
                  response = splotdata$Species_richness,
                  verbose = FALSE,
                  seed = 1,
                  ntree = 5)

  expect_identical(selection$selectedvars, c("bio_6", "bio_12", "bio_15"))
  expect_identical(selection$metric, "RMSE")
  expect_identical(selection$maximize, FALSE)
})


test_that("ffs works for classification with default arguments",{
  data("splotdata")
  splotdata = splotdata |> sf::st_drop_geometry()
  splotdata$Biome = droplevels(splotdata$Biome)
  set.seed(1)
  selection = ffs(predictors = splotdata[,c(6:16)],
                  response = splotdata$Biome,
                  verbose = FALSE,
                  seed = 1,
                  ntree = 5)

  expect_identical(selection$selectedvars, c("bio_4", "bio_6",  "bio_13",
                                             "bio_14", "bio_12", "bio_8", "elev"))
  expect_identical(selection$metric, "Accuracy")
  expect_identical(selection$maximize, TRUE)

})


test_that("ffs works for withinSE = TRUE",{
  data("splotdata")
  splotdata = splotdata |> sf::st_drop_geometry()
  splotdata$Biome = droplevels(splotdata$Biome)
  set.seed(1)
  selection = ffs(predictors = splotdata[,c(6:16)],
                  response = splotdata$Biome,
                  seed = 1,
                  verbose = FALSE,
                  ntree = 5,
                  withinSE = TRUE)

  expect_identical(selection$selectedvars, c("bio_4", "bio_6",  "bio_12",
                                             "bio_14", "bio_8"))


})











## Iris tests that should fail if implemented new


test_that("ffs works with default arguments and the iris dataset",{
  data(iris)
  set.seed(1)
  selection = ffs(predictors = iris[,1:4],
                  response = iris$Species,
                  seed = 1)

  expect_identical(selection$selectedvars, c("Petal.Length", "Petal.Width", "Sepal.Width"))
  expect_equal(selection$selectedvars_perf, c(0.9530141, 0.9544820, 0.9544820))

})



test_that("ffs works with globalVal = TRUE", {

  data(iris)
  set.seed(1)
  selection = ffs(predictors = iris[,1:4],
                  response = iris$Species,
                  seed = 1,
                  globalval = TRUE)

  expect_identical(selection$selectedvars, c("Petal.Length", "Petal.Width", "Sepal.Width"))
  expect_equal(selection$selectedvars_perf, c("Accuracy" = 0.9530792,"Accuracy" = 0.9545455,"Accuracy" = 0.9545455 ), tolerance = 0.005)

})

test_that("ffs works with withinSE = TRUE", {

  data(iris)
  set.seed(1)
  selection = ffs(predictors = iris[,1:4],
                  response = iris$Species,
                  seed = 1,
                  withinSE = TRUE)


  expect_identical(selection$selectedvars, c("Petal.Length", "Petal.Width"))
  expect_equal(selection$selectedvars_perf, c(0.9530141), tolerance = 0.0005)

})


test_that("ffs fails with minvar set to maximum", {

  data(iris)
  set.seed(1)
  expect_error(ffs(predictors = iris[,1:4],
                   response = iris$Species,
                   seed = 1,
                   minVar = 4), regexp = ".*undefined columns selected")



})








