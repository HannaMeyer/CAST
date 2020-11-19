#' Calibrate the AOA based on the relationship between the DI and the prediction error
#' @description
#' @param AOA the result of ?aoa
#' @param model the model used to get the AOA
#' @param calibration Optional. Data frame containing a new set of predictor and response values.
#' @param response character. Name of the response variable used in the model and present in the calibration data frame. Only required when calibration data are given.
#' @param showplot Logical. Show visualization oder not?
#' @details The AOA is the area to which the model can be applied and where the reliability of predictions is expected to be comparable to the cross-validation error of the model.
#' However, it might be desireable to limit predictions to an area with a user-defined error.
#' This function allows for this based on the relationship of the DI an the absolute prediction error derived from cross-validation during model training.
#'
#' The error-DI relationship can also be estimated by providing new external calibration data points.
#' In this case the DI is calculated from these data and the error~DI relationship is derived.
#' If the AOA threshold based on the calibration data is larger than the original AOA threshold, the AOA is updated accordingly.
#' @return rasterStack which contains the original DI and the AOA (whoch might be updated if new test data indicate this option), as well as the expected error.
#'
#' @author
#' Hanna Meyer
#' @seealso \code{\link{aoa}}
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#' library(viridis)
#' library(latticeExtra)
#'
#' # prepare sample data:
#' library(sf)
#' library(raster)
#' library(caret)
#' # prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))[[1:8]]
#' dat <- extract(studyArea,pts,df=TRUE)
#' trainDat <- merge(dat,pts,by.x="ID",by.y="ID")
#'
#' # train a model:
#' variables <- c("DEM","NDRE.Sd","TWI")
#' set.seed(100)
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#'   trainDat$VW,method="rf",importance=TRUE,tuneLength=1,
#'   trControl=trainControl(method="cv",number=5,savePredictions=TRUE))
#'
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea,model)
#'
# '#...or get the area for a user-defined performance
#' AOA_new <- calibrate_aoa(AOA,model)
#' plot(AOA_new$expectedError)
#'
#'
#' ####example with new calibration data:
#' #...to do
#'
#' }
#' @export calibrate_aoa
#' @aliases calibrate_aoa

calibrate_aoa <- function(AOA,model,calibration=NULL,response=NULL,showplot=TRUE){
  ### Get cross-validated predictions from the model:
  preds <- model$pred
  for (i in 1:length(model$bestTune)){
    tunevar <- names(model$bestTune[i])
    preds <- preds[preds[,tunevar]==model$bestTune[,tunevar],]
  }
  preds <- preds[order(preds$rowIndex),]

  ### Estimate the error~DI relationship:
  if(is.null(preds)&is.null(calibration)){
    stop("no cross-predictions can be retrieved from the model. Train with savePredictions=TRUE or provide calibration data")
  }
  if(!is.null(calibration)){
    predictions <- predict(model,calibration)
    DI <- aoa(calibration,model)$DI
    preds <- data.frame("pred"=predictions,
                       "obs"=calibration[,response])


  }  else{

    DI <- attributes(AOA)$TrainDI$DI
  }

  error <- abs(preds$pred-preds$obs)
  #error <- sqrt((preds$pred-preds$obs)^2)

  errormodel <- lm(error~DI)

#  ### Get the threshold for the proficiency value:
#  if(is.null(proficiency)){
#    proficiency <-  predict(errormodel,data.frame("DI"=th_orig))
#  }
#  AOAthres <- (-coefficients(summary(errormodel))[1]+proficiency)/coefficients(summary(errormodel))[2]
#  if (AOAthres<0){
#    AOAthres <- 0
#    message(paste0("Warning: This proficiency cannot be reached. Best proficiency is ",
#                   predict(errormodel,data.frame("DI"=AOAthres))," at a DI-threshold of ",AOAthres))
#  }

  ### Update the AOA
#  AOA$AOA[AOA$DI>AOAthres] <- 0
 # AOA$AOA[AOA$DI<=AOAthres] <- 1
  #attributes(AOA)$aoa_stats$threshold <- AOAthres

  if(!is.null(calibration)){
    AOAthres <- boxplot.stats(DI)$stats[5]
  }else{
    AOAthres <- attributes(AOA)$aoa_stats$threshold
  }

  if(!is.null(calibration)&AOAthres>attributes(AOA)$aoa_stats$threshold){
    attributes(AOA)$aoa_stats$threshold <- AOAthres
    AOA$AOA[AOA$DI>AOAthres] <- 0
    AOA$AOA[AOA$DI<=AOAthres] <- 1
    message("AOA was updated according to the new training data")
  }
  AOA$expectedError <- raster::predict(AOA$DI,errormodel)
  #AOA$expectedError[AOA$AOA==0] <- NA


#  if (inherits(AOA$DI, "Raster")){
#  RMSEv <- predict(AOA$DI,errormodel)
#  RMSEv <- sqrt(mean(values(RMSEv)^2,na.rm=TRUE))
#  }else{
#    RMSEv <- predict(errormodel,AOA$DI)
#    RMSEv <- sqrt(mean(RMSEv^2,na.rm=TRUE))
#  }

  ### Check if the new threshold is based on extrapolation:
 # if(AOAthres>max(DI)){message("Warning: new threshold > maximum DI of the training data. Assumptions about the threshold are made based on extrapolation.")}
  dat <- data.frame(DI,error)

  ### Plot the error~DI relationship and the new threshold
    plt <- ggplot(dat, aes(x = DI, y = error))+geom_point(aes(DI,error,colour="Training data"))+
      scale_color_manual(name="",values = c("Training data" = 'black'))+
      xlab("DI")+ylab("absolute error")+
      xlim(0, max(DI,AOAthres,na.rm=T))+
      stat_smooth(method="lm", fullrange=TRUE)+
      geom_vline(aes(xintercept=attributes(AOA)$aoa_stats$threshold,
                     linetype="AOA threshold"),col="red")+
     # geom_vline(aes(xintercept=AOAthres,linetype="new threshold"),col="red")+
      #scale_linetype_manual(name="",values=c("original threshold"="dashed",
      #                                       "new threshold"="solid"))+
      scale_linetype_manual(name="",values=c("AOA threshold"="dashed"))+
      theme_bw()
    if(showplot){
    suppressMessages(print(plt))
  }
  attributes(AOA)$calibrate_aoa <- plt
#  attributes(AOA)$expectedRMSE <- RMSEv
#  message(paste0("Using the selected proficiency value will lead to an expected RMSE within the AOA of ",
 #                RMSEv))
  return(AOA)
}

