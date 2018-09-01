#' Best subset feature selection
#' @description Evaluate all combinations of predictors during model training
#' @param predictors see \code{\link{train}}
#' @param response see \code{\link{train}}
#' @param method see \code{\link{train}}
#' @param metric see \code{\link{train}}
#' @param maximize see \code{\link{train}}
#' @param trControl see \code{\link{train}}
#' @param tuneLength see \code{\link{train}}
#' @param tuneGrid see \code{\link{train}}
#' @param seed A random number
#' @param verbose Logical. Should information about the progress be printed?
#' @param ... arguments passed to the classification or regression routine
#' (such as randomForest). Errors will occur if values for tuning parameters are
#' passed here.
#' @return A list of class train. Beside of the usual train content
#' the object contains "perf_all" that gives the performance of all model runs.
#' @details bss is an alternative to \code{\link{ffs}} and ideal if the training
#' set is small. Models are iteratively fitted using all different combinations
#' of predictor variables. Hence, 2^X models are calculated. Dont try running bss
#' on very large datasets because the computation time is much higher compared to
#' \code{\link{ffs}}.

#'
#' @note This validation is particulary suitable for
#' leave-one-station-out cross validations where variable selection
#' MUST be based on the performance of the model on the hold out station.
#' A more time efficient alternative is the forward feature selection (\code{\link{ffs}})
#'  (\code{\link{ffs}}).
#' @author Hanna Meyer
#' @seealso \code{\link{train}},\code{\link{ffs}},
#' \code{\link{trainControl}},\code{\link{CreateSpacetimeFolds}}
#' @examples
#' \dontrun{
#' data(iris)
#' bssmodel <- bss(iris[,1:4],iris$Species)
#' bssmodel$finalModel$xNames
#' }
#' @export bss
#' @aliases bss

bss <- function (predictors,
                 response,
                 method = "rf",
                 metric = ifelse(is.factor(response), "Accuracy", "RMSE"),
                 maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                 trControl = caret::trainControl(),
                 tuneLength = 3,
                 tuneGrid = NULL,
                 seed = 100,
                 verbose=TRUE,
                 ...){
  trControl$returnResamp <- "final"
  if(class(response)=="character"){
    response <- factor(response)
    if(metric=="RMSE"){
      metric <- "Accuracy"
      maximize <- TRUE
    }
  }
  se <- function(x){sd(x, na.rm = TRUE)/sqrt(length(na.exclude(x)))}
  n <- length(names(predictors))
  if(maximize) evalfunc <- function(x){max(x,na.rm=T)}
  if(!maximize) evalfunc <- function(x){min(x,na.rm=T)}
  isBetter <- function (actmodelperf,bestmodelperf,maximization=maximize){
    result <- ifelse (!maximization, actmodelperf < bestmodelperf,
                      actmodelperf > bestmodelperf)
    return(result)
  }
  testgrid <- expand.grid(lapply(seq_along(names(predictors)), c, 0))
  testgrid <- testgrid[-which(rowSums(testgrid==0)>=(length(names(predictors))-1)),]
  acc <- 0

  perf_all <- data.frame(matrix(ncol=length(predictors)+3,nrow=nrow(testgrid)))
  names(perf_all) <- c(paste0("var",1:length(predictors)),metric,"SE","nvar")

  for (i in 1:nrow(testgrid)){
    set.seed(seed)
    model <- caret::train(predictors[,unlist(testgrid[i,])],
                   response,method=method,trControl=trControl,
                   tuneLength=tuneLength,
                   tuneGrid=tuneGrid,...)
    actmodelperf <- evalfunc(model$results[,names(model$results)==metric])
    actmodelperfSE <- se(
      sapply(unique(model$resample$Resample),
             FUN=function(x){mean(model$resample[model$resample$Resample==x,
                                                 metric],na.rm=TRUE)}))
    if (i == 1){
      bestmodelperf <- actmodelperf
      bestmodelperfSE <- actmodelperfSE
      bestmodel <- model
    } else{
      if (isBetter(actmodelperf,bestmodelperf,maximization=maximize)){
        bestmodelperf <- actmodelperf
        bestmodelperfSE <- actmodelperfSE
        bestmodel <- model
      }
    }
    acc <- acc+1
    perf_all[acc,1:length(model$finalModel$xNames)] <- model$finalModel$xNames
    perf_all[acc,(length(predictors)+1):ncol(perf_all)] <- c(actmodelperf,actmodelperfSE,length(model$finalModel$xNames))
    if (verbose){
    print(paste0("models that still need to be trained: ",
                 2^n-(n+1) - acc))
    }
  }
  if (maximize){
    selectedvars_perf <- max(bestmodel$results[,metric])
  } else{
    selectedvars_perf <- min(bestmodel$results[,metric])
  }
  bestmodel$selectedvars <- bestmodel$finalModel$xNames
  bestmodel$selectedvars_perf <- selectedvars_perf
  bestmodel$perf_all <- perf_all
  bestmodel$perf_all <- bestmodel$perf_all[!apply(is.na(bestmodel$perf_all), 1, all),]
  bestmodel$perf_all <- bestmodel$perf_all[order(bestmodel$perf_all$nvar),]
  return(bestmodel)
}
