#' Model the relationship between the DI and the prediction error
#' @description Performance metrics are calculated for moving windows of DI values of cross-validated training data
#' @param model the model used to get the AOA
#' @param trainDI the result of \code{\link{trainDI}} or aoa object \code{\link{aoa}}
#' @param multiCV Logical. Re-run model fitting and validation with different CV strategies. See details.
#' @param window.size Numeric. Size of the moving window. See \code{\link{rollapply}}.
#' @param calib Character. Function to model the DI~performance relationship. Currently lm and scam are supported
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer. See ?aoa for further explanation
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @details If multiCV=TRUE the model is re-fitted and validated by length.out new cross-validations where the cross-validation folds are defined by clusters in the predictor space,
#' ranging from three clusters to LOOCV. Hence, a large range of DI values is created during cross-validation.
#' If the AOA threshold based on the calibration data from multiple CV is larger than the original AOA threshold (which is likely if extrapolation situations are created during CV),
#' the AOA threshold changes accordingly. See Meyer and Pebesma (2021) for the full documentation of the methodology.
#' @return A scam or linear model
#' @author
#' Hanna Meyer, Marvin Ludwig
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \doi{10.1111/2041-210X.13650}
#' @seealso \code{\link{aoa}}
#' @example inst/examples/ex_DItoErrormetric.R
#'
#'
#' @export



DItoErrormetric <- function(model, trainDI, multiCV=FALSE,
                            length.out = 10, window.size = 5, calib = "scam",
                            method= "L2", useWeight=TRUE,
                            k = 6, m = 2){


  if(inherits(trainDI,"aoa")){
    trainDI = trainDI$parameters
  }


  # get DIs and Errormetrics OR calculate new ones from multiCV
  if(!multiCV){
    preds_all <- get_preds_all(model, trainDI)
  }
  if(multiCV){
    preds_all <- multiCV(model, length.out, method, useWeight)
  }

  # train model between DI and Errormetric
  error_model = errorModel(preds_all, model, window.size, calib,  k, m)

  # save AOA threshold and raw data
  attr(error_model, "AOA_threshold") <- attr(preds_all, "AOA_threshold")
  class(error_model) <- c("errorModel", class(error_model))
  return(error_model)
}






#' Model expected error between Metric and DI
#' @param preds_all data.frame: pred, obs, DI
#' @param model the model used to get the AOA
#' @param window.size Numeric. Size of the moving window. See \code{\link{rollapply}}.
#' @param calib Character. Function to model the DI~performance relationship. Currently lm and scam are supported
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @return scam or lm
#'

errorModel <- function(preds_all, model, window.size, calib, k, m){

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
  attr(errormodel, "performance") = performance

  return(errormodel)
}


#' MultiCV
#' @description
#' Multiple Cross-Validation with increasing feature space clusteres
#' @param model the model used to get the AOA
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer. See ?aoa for further explanation
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param ... additional parameters to trainDI
#' @returns preds_all
#'
#'

multiCV <- function(model, length.out, method, useWeight,...){

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
    trainDI_new <- trainDI(model_new, method=method, useWeight=useWeight)


    # get cross-validated predictions, order them  and use only those located in the AOA
    preds <- model_new$pred
    preds <- preds[order(preds$rowIndex),c("pred","obs")]
    preds_dat_tmp <- data.frame(preds,"DI"=trainDI_new$trainDI)
    preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$DI <= trainDI_new$threshold,]
    preds_all <- rbind(preds_all,preds_dat_tmp)
  }

  attr(preds_all, "AOA_threshold") <- trainDI_new$threshold
  message(paste0("Note: multiCV=TRUE calculated new AOA threshold of ", round(trainDI_new$threshold, 5),
                 "\nThreshold is stored in the attributes, access with attr(error_model, 'AOA_threshold').",
                 "\nPlease refere to examples and details for further information."))
  return(preds_all)

}


#' Get Preds all
#' @param model, a model
#' @param trainDI, a trainDI
#'

get_preds_all <- function(model, trainDI){

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


  ## add DI from trainDI
  preds_all$DI <- trainDI$trainDI[!is.na(trainDI$trainDI)]
  ## only take predictions from inside the AOA:
  preds_all <-  preds_all[preds_all$DI<=trainDI$threshold,]
  attr(preds_all, "AOA_threshold") <- trainDI$threshold

  return(preds_all)

}

