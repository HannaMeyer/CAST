#' Model the relationship between the DI and LPD and the prediction error
#' @description Performance metrics are calculated for moving windows of DI and LPD values of cross-validated training data
#' @param model the model used to get the AOA
#' @param trainDI the result of \code{\link{trainDI}} or aoa object \code{\link{aoa}}
#' @param multiCV Logical. Re-run model fitting and validation with different CV strategies. See details.
#' @param window.size Numeric. Size of the moving window. See \code{\link{rollapply}}.
#' @param calib Character. Function to model the DI+LPD~performance relationship. Currently lm, scam, lm_exp and scam_exp are supported
#' @param length.out Numeric. Only used if multiCV=TRUE. Number of cross-validation folds. See details.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer. See ?aoa for further explanation
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @details If multiCV=TRUE the model is re-fitted and validated by length.out new cross-validations where the cross-validation folds are defined by clusters in the predictor space,
#' ranging from three clusters to LOOCV. Hence, a large range of DI and DI values is created during cross-validation.
#' If the AOA threshold based on the calibration data from multiple CV is larger than the original AOA threshold (which is likely if extrapolation situations are created during CV),
#' the AOA threshold changes accordingly. See Meyer and Pebesma (2021) for the full documentation of the methodology.
#' @return A scam or linear model
#' @author
#' Hanna Meyer, Marvin Ludwig
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \doi{10.1111/2041-210X.13650}
#' @seealso \code{\link{aoa}}
#' @example inst/examples/ex_DI_LPDtoErrormetric.R
#'
#'
#' @export
#'

DI_LPDtoErrormetric <- function(model, trainDI, multiCV=FALSE,
                            length.out = 10, window.size = 5, calib = "scam",
                            method= "L2", useWeight=TRUE,
                            kDI = 6, mDI = 2, kLPD = 6, mLPD = 2){


  if(inherits(trainDI,"aoa")){
    trainDI = trainDI$parameters
  }


  # get DIs and Errormetrics OR calculate new ones from multiCV
  if(!multiCV){
    preds_all <- get_preds_all_DI_LPD(model, trainDI)
  }
  if(multiCV){
    preds_all <- multiCV_DI_LPD(model, length.out, method, useWeight)
  }

  # train model between DI and Errormetric
  error_model = errorModel_DI_LPD(preds_all, model, window.size, calib,  kDI, mDI, kLPD, mLPD)

  # save AOA threshold and raw data
  attr(error_model, "AOA_threshold") <- attr(preds_all, "AOA_threshold")
  class(error_model) <- c("errorModelDI_LPD", class(error_model))
  return(error_model)
}






#' Model expected error between Metric and DI and LPD
#' @param preds_all data.frame: pred, obs, DI, LPD
#' @param model the model used to get the AOA
#' @param window.size Numeric. Size of the moving window. See \code{\link{rollapply}}.
#' @param calib Character. Function to model the DI~performance relationship. Currently lm and scam are supported
#' @param k Numeric. See mgcv::s
#' @param m Numeric. See mgcv::s
#' @return scam or lm
#'

errorModel_DI_LPD <- function(preds_all, model, window.size, calib, kDI, mDI, kLPD, mLPD){

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

  # performance$metric <- abs(performance$pred - performance$obs)
  # order data according to DI:
  performanceDI <- preds_all[order(preds_all$DI),]
  # calculate performance for moving window:
  performanceDI$metric <- zoo::rollapply(performanceDI[,1:2], window.size,
                                       FUN=function(x){evalfunc(x[,1],x[,2])},
                                       by.column=F,align = "center",fill=NA)
  performanceDI$ll <- data.table::shift(performanceDI$DI,window.size/2)
  performanceDI$ul <- data.table::shift(performanceDI$DI,-round(window.size/2),0)
  performanceDI <- performanceDI[!is.na(performanceDI$metric),]
  performanceDI$ID <- as.numeric(rownames(performanceDI))

  # order data according to LPD:
  performanceLPD <- preds_all[order(preds_all$LPD),]
  # calculate performance for moving window:
  performanceLPD$metric <- zoo::rollapply(performanceLPD[,1:2], window.size,
                                         FUN=function(x){evalfunc(x[,1],x[,2])},
                                         by.column=F,align = "center",fill=NA)
  performanceLPD$ll <- data.table::shift(performanceLPD$LPD,window.size/2)
  performanceLPD$ul <- data.table::shift(performanceLPD$LPD,-round(window.size/2),0)
  performanceLPD <- performanceLPD[!is.na(performanceLPD$metric),]
  performanceLPD$ID <- as.numeric(rownames(performanceLPD))

  performance <- merge(performanceDI, performanceLPD, by = c("ID", "DI", "LPD", "pred", "obs"))
  performance$metric <- (performance$metric.x + performance$metric.y) / 2
  names(performance) <- c("ID", "DI", "LPD", "pred", "obs", "metric.DI", "ll.DI", "ul.DI", "metric.LPD", "ll.LPD", "ul.LPD", "metric")

  # performance <- preds_all[order(preds_all$DI),]
  # # calculate performance for moving window:
  # performance$metric <- zoo::rollapply(performance[,1:2], window.size,
  #                                      FUN=function(x){evalfunc(x[,1],x[,2])},
  #                                      by.column=F,align = "center",fill=NA)
  # performance$ll <- data.table::shift(performance$DI,window.size/2)
  # performance$ul <- data.table::shift(performance$DI,-round(window.size/2),0)
  # performance$LPD <- zoo::rollapply(performance[,4], window.size,
  #                                   FUN=function(x){round(mean(x), digits = 0)},
  #                                   by.column=F,align = "center",fill=NA)
  # performance <- performance[!is.na(performance$metric),]

  ### Estimate Error:
  if(calib=="lm"){
    errormodel <- lm(metric ~ DI+ LPD, data = performance)
  }
  if(calib=="scam"){
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package \"scam\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (model$maximize){ # e.g. accuracy, kappa, r2
      bsDI="mpd"
    }else{
      bsDI="mpi" #e.g. RMSE
    }

    if (model$maximize){ # e.g. accuracy, kappa, r2
      bsLPD="mpi"
    }else{
      bsLPD="mpd" #e.g. RMSE
    }

    errormodel <- scam::scam(metric~s(DI, k=kDI, bs=bsDI, m=mDI)+s(LPD, k=kLPD, bs=bsLPD, m=mLPD),
                             data=performance,
                             family=stats::gaussian(link="identity"))
  }
  if(calib=="lm_exp"){
    errormodel <- lm(metric ~ DI + log(LPD), data = performance)
  }
  if(calib=="scam_exp"){
    if (!requireNamespace("scam", quietly = TRUE)) {
      stop("Package \"scam\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (model$maximize){ # e.g. accuracy, kappa, r2
      bsDI="mpd"
    }else{
      bsDI="mpi" #e.g. RMSE
    }
    errormodel <- scam::scam(metric ~ s(DI, k=kDI, bs=bsDI, m=mDI) + log(LPD),
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

multiCV_DI_LPD <- function(model, length.out, method, useWeight,...){

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
    trainDI_new <- trainDI(model_new, method=method, useWeight=useWeight, LPD = TRUE)


    # get cross-validated predictions, order them  and use only those located in the AOA
    preds <- model_new$pred
    preds <- preds[order(preds$rowIndex),c("pred","obs")]
    preds_dat_tmp <- data.frame(preds,"DI"=trainDI_new$trainDI, "LPD"=trainDI_new$trainLPD)
    preds_dat_tmp <-  preds_dat_tmp[preds_dat_tmp$DI <= trainDI_new$threshold,]
    preds_all <- rbind(preds_all,preds_dat_tmp)
  }

  attr(preds_all, "AOA_threshold") <- trainDI_new$threshold
  attr(preds_all, "avrgLPD") <- trainDI_new$avrgLPD
  message(paste0("Note: multiCV=TRUE calculated new AOA threshold of ", round(trainDI_new$threshold, 5), "and new average LPD of ", trainDI_new$avrgLPD,
                 "\nThreshold is stored in the attributes, access with attr(error_model, 'AOA_threshold').",
                 "\nAverage LPD is stored in the attributes, access with attr(error_model, 'avrgLPD').",
                 "\nPlease refere to examples and details for further information."))
  return(preds_all)

}


#' Get Preds all
#' @param model, a model
#' @param trainDI, a trainDI
#'

get_preds_all_DI_LPD <- function(model, trainDI){

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


  ## add DI from trainDI and LPD from trainLPD
  preds_all$DI <- trainDI$trainDI[!is.na(trainDI$trainDI)]
  preds_all$LPD <- trainDI$trainLPD[!is.na(trainDI$trainLPD)]
  ## only take predictions from inside the AOA:
  preds_all <-  preds_all[preds_all$DI<=trainDI$threshold,]
  # preds_all <-  preds_all[preds_all$LPD>0,]
  attr(preds_all, "AOA_threshold") <- trainDI$threshold
  attr(preds_all, "avrgLPD") <- trainDI$avrgLPD

  return(preds_all)

}

