loaddata <- function() {
  # prepare sample data:
  data(cookfarm)
  dat <- aggregate(cookfarm[,c("VW","Easting","Northing")],by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- sf::st_as_sf(dat,coords=c("Easting","Northing"))
  pts$ID <- 1:nrow(pts)
  set.seed(100)
  pts <- pts[1:30,]
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
  trainDat <- terra::extract(studyArea,pts,na.rm=FALSE)
  trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")

  # train a model:
  set.seed(100)
  variables <- c("DEM","NDRE.Sd","TWI")
  ctrl <- caret::trainControl(method="cv",number=5,savePredictions=T)
  model <- caret::train(trainDat[,which(names(trainDat)%in%variables)],
                 trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
                 trControl=ctrl)


  data <- list(
    studyArea = studyArea,
    trainDat = trainDat,
    variables = variables,
    model = model
  )

  return(data)
}


test_that("AOA works in default: used with raster data and a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 2936)
  # test trainDI
  expect_equal(AOA$parameters$trainDI, c(0.09043580, 0.14046341, 0.16584582, 0.57617177, 0.26840303,
                                         0.14353894, 0.19768329, 0.24022059, 0.06832037, 0.29150668,
                                         0.18471625, 0.57617177, 0.12344463, 0.09043580, 0.14353894,
                                         0.26896008, 0.22713731, 0.24022059, 0.20388725, 0.06832037,
                                         0.23604264, 0.20388725, 0.91513568, 0.09558666, 0.14046341,
                                         0.16214832, 0.37107762, 0.16214832, 0.18471625, 0.12344463))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                "Max.   :4.4485  ", "NA's   :1993  "))
})


test_that("AOA works without a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  ", "NA's   :1993  "))
})

test_that("AOA (including LPD) works with raster data and a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 2936)
  #test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(3, 4, 6, 0, 7,
                                          6, 2, 1, 5, 3,
                                          4, 0, 1, 2, 6,
                                          5, 4, 4, 5, 7,
                                          3, 4, 0, 2, 3,
                                          6, 1, 7, 3, 2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                 "Max.   :4.4485  ", "NA's   :1993  "))
})


test_that("AOA (inluding LPD) works without a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables, LPD = TRUE, maxLPD = 1, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(7, 9, 12, 1, 12,
                                          12, 4, 2, 8, 10,
                                          6, 1, 3,4, 11,
                                          9, 9, 7, 5, 5,
                                          6, 5, 0, 5, 9,
                                          8, 4, 11, 3,2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  ", "NA's   :1993  "))
})


test_that("AOA (including LPD) works in parallel with raster data and a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = F, parallel = TRUE, cores = 2) # limit to 2 cores

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 2936)
  #test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(3, 4, 6, 0, 7,
                                          6, 2, 1, 5, 3,
                                          4, 0, 1, 2, 6,
                                          5, 4, 4, 5, 7,
                                          3, 4, 0, 2, 3,
                                          6, 1, 7, 3, 2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                 "Max.   :4.4485  ", "NA's   :1993  "))
})


test_that("AOA (inluding LPD) works in parallel without a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables, LPD = TRUE, maxLPD = 1, verbose = F, parallel = TRUE, cores = 2) # limit to 2 cores

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(7, 9, 12, 1, 12,
                                          12, 4, 2, 8, 10,
                                          6, 1, 3,4, 11,
                                          9, 9, 7, 5, 5,
                                          6, 5, 0, 5, 9,
                                          8, 4, 11, 3,2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  ", "NA's   :1993  "))
})
