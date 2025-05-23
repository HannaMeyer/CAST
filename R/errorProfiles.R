#' Model and inspect the relationship between the prediction error and measures of dissimilarities and distances
#' @description Performance metrics are calculated for moving windows of dissimilarity values based on cross-validated training data
#' @param model the model used to get the AOA
#' @param trainDI the result of \code{\link{trainDI}} or aoa object \code{\link{aoa}}
#' @param locations Optional. sf object for the training data used in model. Only used if variable=="geodist". Note that they must be in the same order as model$trainingData.
#' @param variable Character. Which dissimilarity or distance measure to use for the error metric. Current options are "DI" or "LPD"
#' @param multiCV Logical. Re-run model fitting and validation with different CV strategies. See details.
#' @param window.size Numeric. Size of the moving window. See \code{\link[zoo]{rollapply}}.
#' @param calib Character. Function to model the DI/LPD~performance relationship. Currently lm and scam are supported
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer. See ?aoa for further explanation
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @details If multiCV=TRUE the model is re-fitted and validated by length.out new cross-validations where the cross-validation folds are defined by clusters in the predictor space,
#' ranging from three clusters to LOOCV. Hence, a large range of dissimilarity values is created during cross-validation.
#' If the AOA threshold based on the calibration data from multiple CV is larger than the original AOA threshold (which is likely if extrapolation situations are created during CV),
#' the AOA threshold changes accordingly. See Meyer and Pebesma (2021) for the full documentation of the methodology.
#' @return A scam, linear model or exponential model
#' @author
#' Hanna Meyer, Marvin Ludwig, Fabian Schumacher
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \doi{10.1111/2041-210X.13650}
#' @seealso \code{\link{aoa}}
#' @examples
#' \dontrun{
#' library(CAST)
#' library(sf)
#' library(terra)
#' library(caret)
#'
#' data(splotdata)
#' predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))
#'
#' model <- caret::train(st_drop_geometry(splotdata)[,6:16], splotdata$Species_richness,
#'    ntree = 10, trControl = trainControl(method = "cv", savePredictions = TRUE))
#'
#' AOA <- aoa(predictors, model, LPD = TRUE, maxLPD = 1)
#'
#' ### DI ~ error
#' errormodel_DI <- errorProfiles(model, AOA, variable = "DI")
#' plot(errormodel_DI)
#' summary(errormodel_DI)
#'
#' expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
#' plot(expected_error_DI)
#'
#' ### LPD ~ error
#' errormodel_LPD <- errorProfiles(model, AOA, variable = "LPD")
#' plot(errormodel_LPD)
#' summary(errormodel_DI)
#'
#' expected_error_LPD = terra::predict(AOA$LPD, errormodel_LPD)
#' plot(expected_error_LPD)
#'
#' ### geodist ~ error
#' errormodel_geodist = errorProfiles(model, locations=splotdata, variable = "geodist")
#' plot(errormodel_geodist)
#' summary(errormodel_DI)
#'
#' dist <- terra::distance(predictors[[1]],vect(splotdata))
#' names(dist) <- "geodist"
#' expected_error_DI <- terra::predict(dist, errormodel_geodist)
#' plot(expected_error_DI)
#'
#'
#' ### with multiCV = TRUE (for DI ~ error)
#' errormodel_DI = errorProfiles(model, AOA, multiCV = TRUE, length.out = 3, variable = "DI")
#' plot(errormodel_DI)
#'
#' expected_error_DI = terra::predict(AOA$DI, errormodel_DI)
#' plot(expected_error_DI)
#'
#' # mask AOA based on new threshold from multiCV
#' mask_aoa = terra::mask(expected_error_DI, AOA$DI > attr(errormodel_DI, 'AOA_threshold'),
#'   maskvalues = 1)
#' plot(mask_aoa)
#' }
#'
#'
#' @export errorProfiles
#' @aliases errorProfiles DItoErrormetric



errorProfiles <- function(model,
                          trainDI=NULL,
                          locations=NULL,
                          variable = "DI",
                          multiCV=FALSE,
                          length.out = 10,
                          window.size = 5,
                          calib = "scam",
                          method= "L2",
                          useWeight=TRUE,
                          k = 6,
                          m = 2){


  if(inherits(trainDI,"aoa")){
    trainDI = trainDI$parameters
  }

  if(!is.null(locations)&variable=="geodist"){
    message("warning: Please ensure that the order of the locations matches to model$trainingData")
  }


  # get DIs and Errormetrics OR calculate new ones from multiCV
  if(!multiCV){
    preds_all <- get_preds_all(model, trainDI, locations, variable)
  }
  if(multiCV){
    preds_all <- multiCV(model, locations, length.out, method, useWeight, variable)
  }

  # train model between DI and Errormetric
  error_model = errorModel(preds_all, model, window.size, calib,  k, m, variable)

  # save AOA threshold and raw data
  attr(error_model, "AOA_threshold") <- attr(preds_all, "AOA_threshold")
  attr(error_model, "variable") <- attr(preds_all, "variable")
  attr(error_model, "metric") <- attr(preds_all, "metric")
  class(error_model) <- c("errorModel", class(error_model))
  return(error_model)
}






# Model expected error between Metric and DI/LPD
errorModel <- function(preds_all, model, window.size, calib, k, m, variable){

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


  # order data according to DI/LPD:
  performance <- preds_all[order(preds_all[,variable]),]
  # calculate performance for moving window:
  performance$metric <- zoo::rollapply(performance[,1:2], window.size,
                                       FUN=function(x){evalfunc(x[,1],x[,2])},
                                       by.column=F,align = "center",fill=NA)
  performance$ll <- data.table::shift(performance[,variable],window.size/2)
  performance$ul <- data.table::shift(performance[,variable],-round(window.size/2),0)
  performance <- performance[!is.na(performance$metric),]

  performance <-  performance[,c(variable,"metric")]
  ### Estimate Error:
  if(calib=="lm"){
    errormodel <- lm(metric ~ ., data = performance)
  }

  if(calib=="scam"){
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package \"scam\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (variable %in% c("DI","geodist")) {
      if (model$maximize){ # e.g. accuracy, kappa, r2
        bs="mpd"
      }else{
        bs="mpi" #e.g. RMSE
      }
      if(variable=="DI"){
        errormodel <- scam::scam(metric~s(DI, k=k, bs=bs, m=m),
                                 data=performance,
                                 family=stats::gaussian(link="identity"))
      }else if (variable=="geodist"){
        errormodel <- scam::scam(metric~s(geodist, k=k, bs=bs, m=m),
                                 data=performance,
                                 family=stats::gaussian(link="identity"))
      }

    } else if (variable == "LPD") {
      if (model$maximize){ # e.g. accuracy, kappa, r2
        bs="mpi"
      }else{
        bs="mpd" #e.g. RMSE
      }
      errormodel <- scam::scam(metric~s(LPD, k=k, bs=bs, m=m),
                               data=performance,
                               family=stats::gaussian(link="identity"))
    }
  }
  if(calib=="exp"){
    if (variable %in% c("DI","geodist")) {
      stop("Exponential model currently only implemented for LPD")
    } else if (variable == "LPD") {
      errormodel <- lm(metric ~ log(LPD), data = performance)
    }
  }

  attr(errormodel, "performance") = performance

  return(errormodel)
}


# MultiCV
multiCV <- function(model, locations, length.out, method, useWeight, variable,...){

  preds_all <- data.frame()
  train_predictors <- model$trainingData[,-which(names(model$trainingData)==".outcome")]
  train_response <- model$trainingData$.outcome

  for (nclst in round(seq(3,nrow(train_predictors), length.out = length.out))){
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
    if (variable == "DI") {
      trainDI_new <- trainDI(model_new, method=method, useWeight=useWeight, verbose =FALSE)
    } else if (variable == "LPD") {
      trainDI_new <- trainDI(model_new, method=method, useWeight=useWeight, LPD = TRUE, verbose =FALSE)
    } else if (variable=="geodist"){
      tmp_gd_new <- CAST::geodist(locations,modeldomain=locations,cvfolds = model$control$indexOut)
      geodist_new <- tmp_gd_new[tmp_gd_new$what=="CV-distances","dist"]

    }


    preds <- model_new$pred
    preds <- preds[order(preds$rowIndex),c("pred","obs")]
    # get cross-validated predictions, order them  and use only those located in the AOA
    if (variable == "DI"){
      preds_dat_tmp <- data.frame(preds,"DI"=trainDI_new$trainDI)
      preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$DI <= trainDI_new$threshold,]
      preds_all <- rbind(preds_all,preds_dat_tmp)
    } else if (variable == "LPD"){
      preds_dat_tmp <- data.frame(preds,"LPD"=trainDI_new$trainLPD)
      preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$LPD > 0,]
      preds_all <- rbind(preds_all,preds_dat_tmp)
    } else if (variable == "geodist"){
      preds_dat_tmp <- data.frame(preds,"geodist"=geodist_new)
      preds_all <- rbind(preds_all,preds_dat_tmp)
      # NO AOA used here
    }
  }
  if(variable=="DI"|variable=="LPD"){
    attr(preds_all, "AOA_threshold") <- trainDI_new$threshold
    message(paste0("Note: multiCV=TRUE calculated new AOA threshold of ", round(trainDI_new$threshold, 5),
                   "\nThreshold is stored in the attributes, access with attr(error_model, 'AOA_threshold').",
                   "\nPlease refere to examples and details for further information."))
  }
  attr(preds_all, "variable") <- variable
  attr(preds_all, "metric") <- model$metric
  return(preds_all)
}


# Get Preds all
get_preds_all <- function(model, trainDI, locations, variable){

  if(is.null(model$pred)){
    stop("no cross-predictions can be retrieved from the model. Train with savePredictions=TRUE or provide calibration data")
  }

  ## extract cv predictions from model
  preds_all <- model$pred
  for (i in 1:length(model$bestTune)){
    tunevar <- names(model$bestTune[i])
    preds_all <- preds_all[preds_all[,tunevar]==model$bestTune[,tunevar],]
  }
  preds_all <- preds_all[order(preds_all$rowIndex),c("pred","obs")]


  if (variable == "DI") {
    ## add DI from trainDI
    preds_all$DI <- trainDI$trainDI[!is.na(trainDI$trainDI)]
    ## only take predictions from inside the AOA:
    preds_all <-  preds_all[preds_all$DI<=trainDI$threshold,]
  } else if (variable == "LPD") {
    ## add LPD from trainLPD
    preds_all$LPD <- trainDI$trainLPD[!is.na(trainDI$trainLPD)]
    ## only take predictions from inside the AOA:
    preds_all <-  preds_all[preds_all$LPD>0,]
  } else if(variable=="geodist"){
    tmp_gd <- CAST::geodist(locations,modeldomain=locations,cvfolds = model$control$indexOut)
    preds_all$geodist <- tmp_gd[tmp_gd$what=="CV-distances","dist"]
  }

  attr(preds_all, "AOA_threshold") <- trainDI$threshold
  attr(preds_all, "variable") <- variable
  attr(preds_all, "metric") <- model$metric
  return(preds_all)
}

