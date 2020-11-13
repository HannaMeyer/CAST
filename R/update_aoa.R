#' Update the AOA based on a user-defined proficiency value
#' @description
#' @param AOA the result of ?aoa
#' @param model the model used to get the AOA
#' @param proficiency numeric. error accepted for the AOA. To first get an overview leave it to NA and select it after inspecting the plot
#' @param calibration Optional. Data frame containing a new set of predictor and response values#
#' @param response character. Name of the response variable used in the model and present in the calibration data frame. Only required when calibration data are given.
#' @param showplot Logical. Show visualization oder not?
#' @details The AOA is the area to which the model can be applied (because it has seen it) and where the reliability of predictions is expected to be comparable to the cross-validation error of the model.
#' However, it might be desireable to limit predictions to an area with a user-defined error.
#' This function allows for this based on the relationship of the DI an the absolute prediction error derived from cross-validation during model training.
#' Based on the selected proficiency, the new DI-threshold is derived from the linear model between the DI and the error.
#' Note that this is only recommended for proficiency values that lead to DI-thresholds smaller than the maximum DI of the training data.
#'
#' The error-DI relationship can also be estimated by providing new
#' @return rasterStack which contains the original DI and the updated AOA. The new threshold is stored in the attributes.
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
#' AOA_new <- update_aoa(AOA,model,proficiency=0.06)
#' spplot(stack(AOA$AOA,AOA_new$AOA),main=c("original","new"))
#'
#'
#' ####example with new calibration data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' calibration <- dat[sample(1:nrow(dat),100),]
#' AOA_new <- update_aoa(AOA,model,proficiency=0.06,calibration=calibration,response="VW")
#'
#' }
#' @export update_aoa
#' @aliases update_aoa

update_aoa <- function(AOA,model,calibration=NULL,response=NULL,proficiency=NULL,showplot=TRUE){
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
    testDI <- aoa(calibration,model)$DI
    preds <- data.frame("pred"=predictions,
                       "obs"=calibration[,response])


  }  else{

    testDI <- attributes(AOA)$TrainDI$DI
  }

  absError <- abs(preds$pred-preds$obs)
  th_orig <- attributes(AOA)$aoa_stats$threshold
  errormodel <- lm(absError~testDI)

  ### Get the threshold for the proficiency value:
  if(is.null(proficiency)){
    proficiency <-  predict(errormodel,data.frame("testDI"=th_orig))
  }
  AOAthres <- (-coefficients(summary(errormodel))[1]+proficiency)/coefficients(summary(errormodel))[2]
  if (AOAthres<0){
    AOAthres <- 0
    message(paste0("Warning: This proficiency cannot be reached. Best proficiency is ",
                   predict(errormodel,data.frame("testDI"=AOAthres))," at a DI-threshold of ",AOAthres))
  }

  ### Update the AOA
  AOA$AOA[AOA$DI>AOAthres] <- 0
  AOA$AOA[AOA$DI<=AOAthres] <- 1
  attributes(AOA)$aoa_stats$threshold <- AOAthres

  ### Check if the new threshold is based on extrapolation:
  if(AOAthres>max(testDI)){message("Warning: new threshold > maximum DI of the training data. Assumptions about the threshold are made based on extrapolation.")}
  dat <- data.frame(testDI,absError)

  ### Plot the error~DI relationship and the new threshold
    plt <- ggplot(dat, aes(x = testDI, y = absError))+geom_point(aes(testDI,absError,colour="Training data"))+
      scale_color_manual(name="",values = c("Training data" = 'black'))+
      xlab("DI")+ylab("absolute error")+
      xlim(0, max(testDI,AOAthres,na.rm=T))+
      stat_smooth(method="lm", fullrange=TRUE)+
      geom_vline(aes(xintercept=th_orig,linetype="original threshold"),col="red")+
      geom_vline(aes(xintercept=AOAthres,linetype="new threshold"),col="red")+
      scale_linetype_manual(name="",values=c("original threshold"="dashed",
                                             "new threshold"="solid"))+
      theme_bw()
    if(showplot){
    suppressMessages(print(plt))
  }
  attributes(AOA)$update_aoa <- plt
  return(AOA)
}

