loaddata <- function() {
  # prepare sample data:
  dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
  dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
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
  model <- caret::train(trainDat[,which(names(trainDat)%in%variables)],
                 trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
                 trControl=caret::trainControl(method="cv",number=5,savePredictions=T))


  data <- list(
    studyArea = studyArea,
    trainDat = trainDat,
    variables = variables,
    model = model
  )

  return(data)
}

test_that("trainDI works in default for a trained model", {
dat <- loaddata()
#...then calculate the DI of the trained model:
DI <- trainDI(model=dat$model, verbose = F)

#test threshold:
expect_equal(as.numeric(round(DI$threshold,5)), 0.38986)
# test trainDI
expect_equal(DI$trainDI, c(0.09043580, 0.14046341, 0.16584582, 0.57617177, 0.26840303,
                           0.14353894, 0.19768329, 0.24022059, 0.06832037, 0.29150668,
                           0.18471625, 0.57617177, 0.12344463, 0.09043580, 0.14353894,
                           0.26896008, 0.22713731, 0.24022059, 0.20388725, 0.06832037,
                           0.23604264, 0.20388725, 0.91513568, 0.09558666, 0.14046341,
                           0.16214832, 0.37107762, 0.16214832, 0.18471625, 0.12344463))
# test summary statistics of the DI
expect_equal(as.numeric(colMeans(DI$train)),
             c(795.4426351,4.0277978,0.2577245))
})

test_that("trainDI (with LPD = TRUE) works in default for a trained model", {
  dat <- loaddata()
  #...then calculate the DI of the trained model:
  DI <- trainDI(model=dat$model, LPD = TRUE, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(DI$threshold,5)), 0.38986)
  #test trainLPD
  expect_identical(DI$trainLPD, as.integer(c(3, 4, 6, 0, 7,
                              6, 2, 1, 5, 3,
                              4, 0, 1, 2, 6,
                              5, 4, 4, 5, 7,
                              3, 4, 0, 2, 3,
                              6, 1, 7, 3, 2)))
  # test summary statistics of the DI
  expect_equal(as.numeric(colMeans(DI$train)),
               c(795.4426351,4.0277978,0.2577245))
})
