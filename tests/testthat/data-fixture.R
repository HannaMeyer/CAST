loaddata <- function() {
  # prepare sample data:
  data("cookfarm")
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