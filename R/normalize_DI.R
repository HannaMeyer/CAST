#' Normalize DI values
#' @description
#' The DI is normalized by the DI threshold to allow for a more straightforward interpretation.
#' A value in the resulting DI larger 1 means that the data are more dissimilar than what has been observed during cross-validation.
#' The returned threshold is adjusted accordingly and is, as a consequence, 1.
#' @param AOA An AOA object
#' @return An object of class \code{aoa}
#' @seealso \code{\link{aoa}}
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(caret)
#'
#' # prepare sample data:
#' data(cookfarm)
#' dat <- aggregate(cookfarm[,c("VW","Easting","Northing")],
#'    by=list(as.character(cookfarm$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' set.seed(100)
#' pts <- pts[1:30,]
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
#' trainDat <- extract(studyArea,pts,na.rm=FALSE)
#' trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
#'
#' # train a model:
#' set.seed(100)
#' variables <- c("DEM","NDRE.Sd","TWI")
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
#' trControl=trainControl(method="cv",number=5,savePredictions=T))
#'
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea, model)
#' plot(AOA)
#' plot(AOA$DI)
#'
#' #... then normalize the DI
#' DI_norm <- normalize_DI(AOA)
#' plot(DI_norm)
#' plot(DI_norm$DI)
#'
#' }
#' @export normalize_DI
#' @aliases normalize_DI


normalize_DI <- function(AOA) {
  AOA$DI <- AOA$DI/AOA$parameters$threshold
  AOA$parameters$trainDI <- AOA$parameters$trainDI/AOA$parameters$threshold
  AOA$parameters$threshold <- AOA$parameters$threshold/AOA$parameters$threshold
  return(AOA)
}

