
test_that("errorProfiles works in default settings", {
  data(splotdata)
  splotdata <- sf::st_drop_geometry(splotdata)
  predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

  set.seed(100)
  model <- caret::train(splotdata[,6:16], splotdata$Species_richness, ntree = 10,
                        trControl = caret::trainControl(method = "cv", savePredictions = TRUE))

  AOA <- CAST::aoa(predictors, model)

  # DI ~ error
  errormodel_DI <- CAST::errorProfiles(model, AOA, variable = "DI")

  expected_error_DI = terra::predict(AOA$DI, errormodel_DI)



  #test model fit:
  expect_equal(round(as.numeric(summary(errormodel_DI$fitted.values)),2),
               c(14.25, 14.34, 15.21, 17.23, 18.70, 27.46))
  # test model predictions
  expect_equal(as.vector( summary(terra::values(expected_error_DI))),
               c("Min.   :14.26  ", "1st Qu.:27.46  ", "Median :27.46  ",
                 "Mean   :26.81  ", "3rd Qu.:27.46  ","Max.   :27.47  ",
                 "NA's   :17678  "))
})


test_that("errorProfiles works in with LPD", {
  data(splotdata)
  splotdata <- sf::st_drop_geometry(splotdata)
  predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

  set.seed(100)
  model <- caret::train(splotdata[,6:16], splotdata$Species_richness, ntree = 10,
                        trControl = caret::trainControl(method = "cv", savePredictions = TRUE))

  AOA <- CAST::aoa(predictors, model, LPD = TRUE, maxLPD = 1)
  errormodel_LPD <- CAST::errorProfiles(model, AOA, variable = "LPD")
  expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)


  #test model fit:
  expect_equal(round(as.numeric(summary(errormodel_LPD$fitted.values)),2),
               c(16.36, 16.36, 16.36, 16.36, 16.36, 16.36))
  # test model predictions
  expect_equal(as.vector(summary(terra::values(expected_error_LPD))),
               c("Min.   :16.36  ", "1st Qu.:16.36  ", "Median :16.36  ",
                 "Mean   :16.36  ", "3rd Qu.:16.36  ",
                 "Max.   :16.36  ", "NA's   :17678  "))

})



test_that("errorProfiles works for multiCV", {
  data(splotdata)
  splotdata <- sf::st_drop_geometry(splotdata)
  predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

  set.seed(100)
  model <- caret::train(splotdata[,6:16], splotdata$Species_richness, ntree = 10,
                        trControl = caret::trainControl(method = "cv", savePredictions = TRUE))

  AOA <- CAST::aoa(predictors, model)
  # with multiCV = TRUE (for DI ~ error)
  set.seed(100)
  errormodel_DI = suppressWarnings(errorProfiles(model, AOA, multiCV = TRUE, length.out = 3))
  expected_error_DI = terra::predict(AOA$DI, errormodel_DI)

  #test model fit:
  expect_equal(round(as.numeric(summary(errormodel_DI$fitted.values)),2),
               c(12.53, 17.21, 26.80, 26.19, 35.28, 35.30))
  # test model predictions
  expect_equal(as.vector( summary(terra::values(expected_error_DI))),
               c("Min.   :13.11  ", "1st Qu.:32.58  ", "Median :35.05  ",
                 "Mean   :32.54  ", "3rd Qu.:35.30  ",
                 "Max.   :35.30  ", "NA's   :17678  "))


})
