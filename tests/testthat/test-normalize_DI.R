
test_that("normalize_DI works with default arguments",{
  skip_if_not_installed("randomForest")
  data(cookfarm)
  dat <- aggregate(cookfarm[,c("VW","Easting","Northing")],
                   by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- sf::st_as_sf(dat,coords=c("Easting","Northing"))
  pts$ID <- 1:nrow(pts)
  pts <- pts[1:30,]
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
  trainDat <- terra::extract(studyArea,pts,na.rm=FALSE)
  trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
  set.seed(100)
  variables <- c("DEM","NDRE.Sd","TWI")
  model <- caret::train(trainDat[,which(names(trainDat)%in%variables)],
                 trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
                 trControl=trainControl(method="cv",number=5,savePredictions=T))
  AOA <- aoa(studyArea, model)
  DI_norm <- normalize_DI(AOA)
  expect_equal(as.integer(DI_norm$parameters$threshold), 1)
  expect_equal(as.vector(summary(terra::values(DI_norm$DI)))[1:6],
               c("Min.   : 0.0000  ", "1st Qu.: 0.3410  ", "Median : 0.5263  ",
                 "Mean   : 0.7330  ", "3rd Qu.: 0.9785  ",
                 "Max.   :11.4105  "))

})
