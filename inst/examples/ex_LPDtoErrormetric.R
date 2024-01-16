\dontrun{

  library(CAST)
  library(sf)
  library(terra)
  library(caret)


  data(splotdata)
  splotdata <- st_drop_geometry(splotdata)
  predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

  model <- caret::train(splotdata[,6:16], splotdata$Species_richness, ntree = 10,
                        trControl = trainControl(method = "cv", savePredictions = TRUE))

  AOA <- aoa(predictors, model, LPD = TRUE, maxLPD = 1)

  errormodel <- LPDtoErrormetric(model, AOA)
  plot(errormodel)

  expected_error = terra::predict(AOA$LPD, errormodel)
  plot(expected_error)


  # with multiCV = TRUE
  errormodel = LPDtoErrormetric(model, AOA, multiCV = TRUE, length.out = 3)
  plot(errormodel)

  expected_error = terra::predict(AOA$LPD, errormodel)
  plot(expected_error)

  # mask AOA based on new threshold from multiCV
  mask_aoa = terra::mask(expected_error, AOA$DI > attr(errormodel, 'AOA_threshold'), maskvalues = 1)
  plot(mask_aoa)




}
