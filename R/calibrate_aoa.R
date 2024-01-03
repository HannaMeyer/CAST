#' Calibrate the AOA based on the relationship between the DI and the prediction error
#' @description Performance metrics are calculated for moving windows of DI values of cross-validated training data
#' @param AOA the result of \code{\link{aoa}}
#' @param model the model used to get the AOA
#' @param window.size Numeric. Size of the moving window. See \code{\link{rollapply}}.
#' @param calib Character. Function to model the DI~performance relationship. Currently lm and scam are supported
#' @param multiCV Logical. Re-run model fitting and validation with different CV strategies. See details.
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param maskAOA Logical. Should areas outside the AOA set to NA?
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer. See ?aoa for further explanation
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @param showPlot Logical.
#' @details If multiCV=TRUE the model is re-fitted and validated by length.out new cross-validations where the cross-validation folds are defined by clusters in the predictor space,
#' ranging from three clusters to LOOCV. Hence, a large range of DI values is created during cross-validation.
#' If the AOA threshold based on the calibration data from multiple CV is larger than the original AOA threshold (which is likely if extrapolation situations are created during CV),
#' the AOA is updated accordingly. See Meyer and Pebesma (2021) for the full documentation of the methodology.
#' @return A list of length 2 with the elements "AOA": SpatRaster or stars object which contains the original DI and the AOA (which might be updated if new test data indicate this option), as well as the expected performance based on the relationship.
#' Data used for calibration are stored in the attributes. The second element is a plot showing the relationship.
#' @author
#' Hanna Meyer
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \doi{10.1111/2041-210X.13650}
#' @seealso \code{\link{aoa}}
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(caret)
#' library(viridis)
#' library(latticeExtra)
#'
#' #' # prepare sample data:
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
#' dat <- extract(studyArea,pts,na.rm=TRUE)
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
# '# and get the expected performance on a pixel-level:
#' AOA_new <- calibrate_aoa(AOA,model)
#' plot(AOA_new$AOA$expected_RMSE)
#' }
#' @export calibrate_aoa
#' @aliases calibrate_aoa

calibrate_aoa <- function(AOA,model, window.size=5, calib="scam",multiCV=FALSE,
                          length.out = 10, maskAOA=TRUE, method= "L2", useWeight=TRUE,
                          showPlot=TRUE,k=6,m=2){


  message("Note: calibrate_aoa is deprecated and will be removed soon. Please use and refere to DItoErrormetric instead.")

  as_stars <- FALSE
#  as_raster <- FALSE

  if (inherits(AOA$AOA, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    attr <- attributes(AOA)[c("aoa_stats","TrainDI")]
    AOA$AOA <- methods::as(AOA$AOA, "SpatRaster")
    AOA$DI <- methods::as(AOA$DI, "SpatRaster")
    attributes(AOA)<- c(attributes(AOA),attr)
    as_stars <- TRUE
  }

  if (inherits(AOA$AOA, "Raster")) {
    #if (!requireNamespace("raster", quietly = TRUE))
    #  stop("package raster required: install that first")
    message("Raster will soon not longer be supported. Use terra or stars instead")
    attr <- attributes(AOA)[c("aoa_stats","TrainDI")]
    AOA$AOA <- methods::as(AOA$AOA, "SpatRaster")
    AOA$DI <- methods::as(AOA$DI, "SpatRaster")
    attributes(AOA)<- c(attributes(AOA),attr)
#    as_raster <- TRUE
  }

  if(multiCV){
    preds_all <- data.frame()
    train_predictors <- model$trainingData[,-which(names(model$trainingData)==".outcome")]
    train_response <- model$trainingData$.outcome

    for (nclst in round(seq(3,nrow(train_predictors),length.out = length.out))){
      # define clusters in predictor space used for CV:
      clstrID <- tryCatch({stats::kmeans(train_predictors,nclst)$cluster},
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
      mcall$method <- model$method
      mcall$metric <- model$metric
      mcall$cl <- NULL # fix option for parallel later

      # retrain model and calculate AOA
      model_new <- do.call(caret::train,mcall)
      AOA_new <- aoa(train_predictors,model_new,method=method, useWeight=useWeight)

      # legacy change (very dirty, change this as soon as possible)
      #AOA_new <- AOA_new$AOA


      # get cross-validated predictions, order them  and use only those located in the AOA
      preds <- model_new$pred
      preds <- preds[order(preds$rowIndex),c("pred","obs")]
      preds_dat_tmp <- data.frame(preds,"DI"=AOA_new$parameters$trainDI)
      preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$DI<=AOA_new$parameters$threshold,]
      preds_all <- rbind(preds_all,preds_dat_tmp)


    }
    AOA$parameters$threshold <- max(preds_all$DI)
    AOA$parameters$trainDI <- preds_all$DI
  }


  ### Get cross-validated predictions from the model:
  if(!multiCV){
    # Get cross-validated predictions:
    if(is.null(model$pred)){
      stop("CV predictions cannot be derived from the model. re-train using savePredictions, see ?caret::trainControl")
    }
    preds_all <- model$pred
    for (i in 1:length(model$bestTune)){
      tunevar <- names(model$bestTune[i])
      preds_all <- preds_all[preds_all[,tunevar]==model$bestTune[,tunevar],]
    }

    preds_all <- preds_all[order(preds_all$rowIndex),c("pred","obs")]
    preds_all$DI <- AOA$parameters$trainDI[!is.na(AOA$parameters$trainDI)]

    ## only take predictions from inside the AOA:
    preds_all <-  preds_all[preds_all$DI<=AOA$parameters$threshold,]
  }


  ### Estimate the error~DI relationship:
  if(is.null(preds_all)){
    stop("no cross-predictions can be retrieved from the model. Train with savePredictions=TRUE or provide calibration data")
  }
  ## use performance metric from the model:
  rmse <- function(pred,obs){sqrt( mean((pred - obs)^2, na.rm = TRUE) )}
  rsquared <-  function(pred,obs){summary(lm(pred~obs))$r.squared}
  mae <- function(pred,obs){MAE(pred,obs)}
  kappa <- function(pred,obs){
    pred <- factor(pred)
    obs <- factor(obs)
    lev <- unique(c(levels(pred), levels(obs)))
    pred <- factor(pred, levels = lev)
    obs <- factor(obs, levels = lev)
    result <- tryCatch( confusionMatrix(pred, obs)$overall["Kappa"], error = function(e)e)
    if(inherits(result, "error")){result <- 0} # 0 not right value!!! adjust!!!
    return(unname(result))
  }

  accuracy <- function(pred,obs){
    pred <- factor(pred)
    obs <- factor(obs)
    lev <- unique(c(levels(pred), levels(obs)))
    pred <- factor(pred, levels = lev)
    obs <- factor(obs, levels = lev)
    result <- tryCatch(confusionMatrix(pred, obs)$overall["Accuracy"], error = function(e)e)
    if(inherits(result, "error")){result <- 0}
    return(unname(result))
  }
  if(!tolower(model$metric)%in%c("rmse","rsquared","mae","kappa","accuracy")){
    message("Model metric not yet included in this function")
    stop()
  }

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

  performance <- performance[!is.na(performance$metric),]

  ### Update AOA:
  if (multiCV){
    if(inherits(AOA$AOA,"SpatRaster")){
      AOA$AOA <- terra::setValues(AOA$AOA, 0)
    }else{
      AOA$AOA <- 0
    }
    AOA$AOA[AOA$DI<=max(performance$DI,na.rm=T)] <- 1
    if(inherits(AOA$AOA,"SpatRaster")){
      AOA$AOA <- terra::mask(AOA$AOA,AOA$DI)
    }else{
      AOA$AOA[is.na(AOA$DI)] <- NA
    }
  }

  ### Estimate Error:
  if(calib=="lm"){
    errormodel <- lm(metric ~ DI, data = performance)
  }
  if(calib=="scam"){
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package \"scam\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (model$maximize){ # e.g. accuracy, kappa, r2
      bs="mpd"
    }else{
      bs="mpi" #e.g. RMSE
    }

    errormodel <- scam::scam(metric~s(DI, k=k, bs=bs, m=m),
                             data=performance,
                             family=stats::gaussian(link="identity"))
  }


  attributes(AOA)$calib$model <- errormodel

  DI_pred <- AOA$DI
  attr <- attributes(AOA)[c("aoa_stats","TrainDI","calib")]
  # predict and make sure it's not going beyond min observed values
  if(inherits(DI_pred,"SpatRaster")){
    terra::values(DI_pred)[terra::values(AOA$DI)<min(performance$DI,na.rm=TRUE)] <- min(performance$DI,na.rm=TRUE)
    AOA$expectedError <- terra::predict(DI_pred,errormodel)
  }else{
    DI_pred[AOA$DI<min(performance$DI,na.rm=TRUE)] <- min(performance$DI,na.rm=TRUE)
    AOA$expectedError <- predict(errormodel,data.frame("DI"=DI_pred))
  }


  if(maskAOA){
    if(inherits(AOA$expectedError,"SpatRaster")){
      AOA$expectedError <-  terra::mask(AOA$expectedError,AOA$AOA,maskvalue=0)
    }else{
      AOA$expectedError[AOA$AOA==0] <- NA
    }
  }

  names(performance)[which(names(performance)=="metric")] <- model$metric
  attr$calib$group_stats <- performance
  attributes(AOA)<- c(attributes(AOA),attr)
  ### Plot result:

  # if(showPlot){

  # loc <- "topleft"
  # if(model$maximize){
  #   loc <- "topright"
  # }
  #
  # plot(attr$calib$group_stats$DI,attr$calib$group_stats[,model$metric],xlab="DI",
  #      ylab=model$metric)
  # graphics::legend(loc,lty=c(NA,2),lwd=c(NA,1),pch=c(1,NA),col=c("black","black"),
  #        legend=c("CV","model"),bty="n")
  # graphics::lines(seq(0,max(attr$calib$group_stats$DI, na.rm=TRUE),max(attr$calib$group_stats$DI, na.rm=TRUE)/100),
  #       predict(attributes(AOA)$calib$model,
  #               data.frame("DI"=seq(0, max(attr$calib$group_stats$DI,na.rm=TRUE),
  #                                  max(attr$calib$group_stats$DI, na.rm=TRUE)/100))),lwd=1,lty=2,col="black")


  p <- lattice::xyplot(attr$calib$group_stats[,model$metric]~attr$calib$group_stats$DI,xlab="DI",
                       ylab=model$metric,col="black",
                       key=list(columns=2,
                                text=list(lab=c("cross-validation","model")),
                                points=list(pch=c(1,NA), col="black"),
                                lines=list(lty=c(0,2), lwd=2, col="black")),panel = function(x, y, ...) {
                                  lattice::panel.xyplot(x, y, ...)
                                  lattice::llines(x, predict(attr$calib$model), col="black", lwd=2, lty=2)
                                })

  if(showPlot){
    print(p)
  }


  if (as_stars){
    AOA$AOA <- split(stars::st_as_stars(AOA$AOA), "band")
    AOA$DI <- split(stars::st_as_stars(AOA$DI), "band")
    AOA$expectedError <- split(stars::st_as_stars(AOA$expectedError), "band")
    attributes(AOA$AOA)<- c(attributes(AOA$AOA),attr)
  }

#  if(as_raster){
#    AOA$AOA <- methods::as(AOA$AOA, "Raster")
#    AOA$DI <- methods::as(AOA$DI, "Raster")
#    AOA$expectedError <- methods::as(AOA$expectedError, "Raster")
#    attributes(AOA$AOA)<- c(attributes(AOA$AOA),attr)
#  }
  names(AOA)[names(AOA)=="expectedError"] <- paste0("expected_",model$metric)
  #return(AOA)

  return(list(AOA = AOA,
              plot = p))

}


