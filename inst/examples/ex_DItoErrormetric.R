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

  # DI ~ error
  errormodel_DI <- DItoErrormetric(model, AOA, what = "DI")
  plot(errormodel_DI)

  expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
  plot(expected_error_DI)

  # LPD ~ error
  errormodel_LPD <- DItoErrormetric(model, AOA, what = "LPD")
  plot(errormodel_LPD)

  expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)
  plot(expected_error_LPD)


  # with multiCV = TRUE (for DI ~ error)
  errormodel_DI = DItoErrormetric(model, AOA, multiCV = TRUE, length.out = 3, what = "DI")
  plot(errormodel_DI)

  expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
  plot(expected_error_DI)

  # mask AOA based on new threshold from multiCV
  mask_aoa = terra::mask(expected_error_DI, AOA$DI > attr(errormodel_DI, 'AOA_threshold'), maskvalues = 1)
  plot(mask_aoa)




}

