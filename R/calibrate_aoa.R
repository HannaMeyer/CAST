#' Calibrate the AOA based on the relationship between the DI and the prediction error
#' @description Performance metrics are calculated for moving windows of DI values of cross-validated training data
#' @param AOA the result of ?aoa
#' @param model the model used to get the AOA
#' @param window.size Numeric. Size of the moving window. See ?zoo::rollapply
#' @param calib Character. Function to model the DI~performance relationship. Currently lm and scam are supported
#' @param multiCV Logical. Re-run model fitting and validation with different CV strategies. See details.
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param maskAOA Logical. Should areas ouside the AOA set to NA?
#' @param showPlot Logical
#' @details The AOA is the area to which the model can be applied and where the reliability of predictions is expected to be comparable to the cross-validation error of the model.
#' However, it might be desireable to limit predictions to an area with a user-defined error.
#' This function allows for this based on the relationship of the DI and the prediction error derived from cross-validation during model training.
#' If multiCV=TRUE the model is re-fitted and validated by length.out new cross-validations where the cross-validation folds are defined by clusters in the predictor space,
#' ranging from three clusters to LOOCV.
#' If the AOA threshold based on the calibration data is larger than the original AOA threshold, the AOA is updated accordingly.
#' @return rasterStack which contains the original DI and the AOA (which might be updated if new test data indicate this option), as well as the expected error based on the relationship. Data used for calibration are stored in the attributes.
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
#' }
#' @export calibrate_aoa
#' @aliases calibrate_aoa

calibrate_aoa <- function(AOA,model, window.size=5, calib="scam",multiCV=FALSE,
                          length.out = 10, maskAOA=TRUE, showPlot=TRUE){

  if(multiCV){
    preds_all <- data.frame()
    train_predictors <- model$trainingData[,-which(names(model$trainingData)==".outcome")]
    train_response <- model$trainingData$.outcome

    for (nclst in round(seq(3,nrow(train_predictors),length.out = length.out))){
      # define clusters in predictor space used for CV:
      clstrID <- tryCatch({kmeans(train_predictors,nclst)$cluster},
                          error=function(e)e)
      if(inherits(clstrID,"error")){next}
      clstrID <- clstrID
      folds <- CreateSpacetimeFolds(data.frame("clstrID"=clstrID), spacevar="clstrID",k=nclst)

      # update model call with new CV strategy:
      mcall <- as.list(model$call)
      mcall <- mcall[-which(names(mcall)%in%c("form","data","","x","y","","trControl"))]
      mcall$x <- quote(train_predictors)
      mcall$y <- quote(train_response)
      mcall$trControl <- trainControl(method="cv",index=folds$index,savePredictions = TRUE)
      mcall$tuneGrid <- model$bestTune

      # retrain model and calculate AOA
      model_new <- do.call(caret::train,mcall)
      AOA_new <- aoa(train_predictors,model_new)

      # get cross-validated predictions, order them  and use only those located in the AOA
      preds <- model_new$pred
      preds <- preds[order(preds$rowIndex),c("pred","obs")]
      preds_dat_tmp <- data.frame(preds,"DI"=attributes(AOA_new)$TrainDI)
      preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$DI<=attributes(AOA_new)$aoa_stats$threshold,]
      preds_all <- rbind(preds_all,preds_dat_tmp)
    }
  }


  ### Get cross-validated predictions from the model:
  if(!multiCV){
    # Get cross-validated predictions:
    preds_all <- model$pred
    for (i in 1:length(model$bestTune)){
      tunevar <- names(model$bestTune[i])
      preds_all <- preds_all[preds_all[,tunevar]==model$bestTune[,tunevar],]
    }

    preds_all <- preds_all[order(preds_all$rowIndex),c("pred","obs")]
    preds_all$DI <- attributes(AOA)$TrainDI

    ## only take predictions from inside the AOA:
    preds_all <-  preds_all[preds_all$DI<=attributes(AOA)$aoa_stats$threshold,]
  }


  ### Estimate the error~DI relationship:
  if(is.null(preds_all)){
    stop("no cross-predictions can be retrieved from the model. Train with savePredictions=TRUE or provide calibration data")
  }
  ## use performance metric from the model:
  evalfunc <- function(pred,obs){
    eval(parse(text=paste0(tolower(model$metric),"(pred,obs)")))
  }

  # order data according to DI:
  performance <- preds_all[order(preds_all$DI),]

  # calculate performance for moving window:
  performance$metric <- zoo::rollapply(performance[,1:2], window.size,
                                FUN=function(x){evalfunc(x[,1],x[,2])},
                    by.column=F,align = "center",fill=NA)

  performance$ll <- data.table::shift(performance$DI,window.size/2)
  performance$ul <- data.table::shift(performance$DI,-round(window.size/2),0)
  #####################################

  #errors <-c()
  #for (i in unique(preds_all$group)){
  #  errors <- c(errors, evalfunc(preds_all$pred[preds_all$group==i],
  #                         preds_all$obs[preds_all$group==i]))

  #}

 # grouped_res <- data.frame("group"=unique(preds_all$group),"performance"=errors)
 # min_max <- unlist(strsplit(gsub("(?![,.])[[:punct:]]", "", as.character(grouped_res$group), perl=TRUE), ","))
 # recl <- data.frame("min"=as.numeric(min_max[seq(1, length(min_max), by=2)]),
   #                  "max"=as.numeric(min_max[seq(2, length(min_max), by=2)]),
   #                  "perf"= grouped_res$performance)
 # recl$DI <- rowMeans(recl[,c("min","max")])


  ### Update AOA:
  if (multiCV){
  AOA$AOA <- 0
  AOA$AOA[AOA$DI<=max(performance$DI,na.rm=T)] <- 1
  if(class(AOA$AOA)=="Raster"){
  AOA$AOA <- raster::mask(AOA$AOA,AOA$DI)
  }else{
    AOA$AOA[is.na(AOA$DI)] <- NA
  }
}

  ### Estimate Error:
if(calib=="lm"){
  errormodel <- lm(metric ~ DI, data = performance)
}
  if(calib=="scam"){
  errormodel <- scam::scam(metric~s(DI, k=5, bs="mpi", m=2),
                    data=performance,
                    family=gaussian(link="identity"))
}


  attributes(AOA)$calib$model <- errormodel

  DI_pred <- AOA$DI

  # predict and make sure it's not going beyond min observed values
  if(class(DI_pred)=="Raster"){
  raster::values(DI_pred)[raster::values(AOA$DI)<min(performance$DI,na.rm=TRUE)] <- min(performance$DI,na.rm=TRUE)
  AOA$expectedError <- raster::predict(DI_pred,errormodel)
  }else{
    DI_pred[AOA$DI<min(performance$DI,na.rm=TRUE)] <- min(performance$DI,na.rm=TRUE)
    AOA$expectedError <- predict(errormodel,data.frame("DI"=DI_pred))
}


  if(maskAOA){
    if(class( AOA$expectedError)=="Raster"){
      AOA$expectedError <-  raster::mask(AOA$expectedError,AOA$AOA,maskvalue=0)
    }else{
      AOA$expectedError[AOA$AOA==0] <- NA
    }
  }

  ### Plot result:

  plot(performance$DI,performance$metric,
       xlab="DI",ylab=model$metric,type="l")
  lines(seq(0,max(performance$DI),0.01),
        predict(errormodel,data.frame("DI"=seq(0,max(performance$DI),0.01))),
        col="red")

  names(performance)[which(names(performance)=="metric")] <- model$metric
  attributes(AOA)$calib$group_stats <- performance
  return(AOA)
}

