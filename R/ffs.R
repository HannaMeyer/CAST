#' Forward feature selection
#' @description A simple forward feature selection algorithm
#' @param predictors see \code{\link{train}}
#' @param response see \code{\link{train}}
#' @param method see \code{\link{train}}
#' @param metric see \code{\link{train}}
#' @param maximize see \code{\link{train}}
#' @param globalval Logical. Should models be evaluated based on 'global' performance? See \code{\link{global_validation}}
#' @param withinSE Logical Models are only selected if they are better than the
#' currently best models Standard error
#' @param minVar Numeric. Number of variables to combine for the first selection.
#' See Details.
#' @param trControl see \code{\link{train}}
#' @param tuneLength see \code{\link{train}}
#' @param tuneGrid see \code{\link{train}}
#' @param seed A random number used for model training
#' @param verbose Logical. Should information about the progress be printed?
#' @param ... arguments passed to the classification or regression routine
#' (such as randomForest).
#' @return A list of class train. Beside of the usual train content
#' the object contains the vector "selectedvars" and "selectedvars_perf"
#' that give the order of the best variables selected as well as their corresponding
#' performance (starting from the first two variables). It also contains "perf_all"
#' that gives the performance of all model runs.
#' @details Models with two predictors are first trained using all possible
#' pairs of predictor variables. The best model of these initial models is kept.
#' On the basis of this best model the predictor variables are iteratively
#' increased and each of the remaining variables is tested for its improvement
#' of the currently best model. The process stops if none of the remaining
#' variables increases the model performance when added to the current best model.
#'
#' The internal cross validation can be run in parallel. See information
#' on parallel processing of carets train functions for details.
#'
#' Using withinSE will favour models with less variables and
#' probably shorten the calculation time
#'
#' Per Default, the ffs starts with all possible 2-pair combinations.
#' minVar allows to start the selection with more than 2 variables, e.g.
#' minVar=3 starts the ffs testing all combinations of 3 (instead of 2) variables
#' first and then increasing the number. This is important for e.g. neural networks
#' that often cannot make sense of only two variables. It is also relevant if
#' it is assumed that the optimal variables can only be found if more than 2
#' are considered at the same time.
#'
#' @note This variable selection is particularly suitable for spatial
#' cross validations where variable selection
#' MUST be based on the performance of the model for predicting new spatial units.
#' See Meyer et al. (2018) and Meyer et al. (2019) for further details.
#'
#' @author Hanna Meyer
#' @seealso \code{\link{train}},\code{\link{bss}},
#' \code{\link{trainControl}},\code{\link{CreateSpacetimeFolds}},\code{\link{nndm}}
#' @references
#' \itemize{
#' \item Gasch, C.K., Hengl, T., Gräler, B., Meyer, H., Magney, T., Brown, D.J. (2015): Spatio-temporal interpolation of soil water, temperature, and electrical conductivity in 3D+T: the Cook Agronomy Farm data set. Spatial Statistics 14: 70-90.
#' \item Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018): Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software 101: 1-9.  \doi{10.1016/j.envsoft.2017.12.001}
#' \item Meyer, H., Reudenbach, C., Wöllauer, S., Nauss, T. (2019): Importance of spatial predictor variable selection in machine learning applications - Moving from data reproduction to spatial prediction. Ecological Modelling. 411, 108815. \doi{10.1016/j.ecolmodel.2019.108815}.
#' \item Ludwig, M., Moreno-Martinez, A., Hölzel, N., Pebesma, E., Meyer, H. (2023): Assessing and improving the transferability of current global spatial prediction models. Global Ecology and Biogeography. \doi{10.1111/geb.13635}.
#' }
#' @examples
#' \dontrun{
#' data(iris)
#' ffsmodel <- ffs(iris[,1:4],iris$Species)
#' ffsmodel$selectedvars
#' ffsmodel$selectedvars_perf
#'}
#'
#' # or perform model with target-oriented validation (LLO CV)
#' #the example is described in Gasch et al. (2015). The ffs approach for this dataset is described in
#' #Meyer et al. (2018). Due to high computation time needed, only a small and thus not robust example
#' #is shown here.
#'
#' \dontrun{
#' #run the model on three cores:
#' library(doParallel)
#' library(lubridate)
#' cl <- makeCluster(3)
#' registerDoParallel(cl)
#'
#' #load and prepare dataset:
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' trainDat <- dat[dat$altitude==-0.3&year(dat$Date)==2012&week(dat$Date)%in%c(13:14),]
#'
#' #visualize dataset:
#' ggplot(data = trainDat, aes(x=Date, y=VW)) + geom_line(aes(colour=SOURCEID))
#'
#' #create folds for Leave Location Out Cross Validation:
#' set.seed(10)
#' indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID",k=3)
#' ctrl <- trainControl(method="cv",index = indices$index)
#'
#' #define potential predictors:
#' predictors <- c("DEM","TWI","BLD","Precip_cum","cday","MaxT_wrcc",
#' "Precip_wrcc","NDRE.M","Bt","MinT_wrcc","Northing","Easting")
#'
#' #run ffs model with Leave Location out CV
#' set.seed(10)
#' ffsmodel <- ffs(trainDat[,predictors],trainDat$VW,method="rf",
#' tuneLength=1,trControl=ctrl)
#' ffsmodel
#' plot(ffsmodel)
#' #or only selected variables:
#' plot(ffsmodel,plotType="selected")
#'
#' #compare to model without ffs:
#' model <- train(trainDat[,predictors],trainDat$VW,method="rf",
#' tuneLength=1, trControl=ctrl)
#' model
#' stopCluster(cl)
#'}
#' @export ffs
#' @aliases ffs

ffs <- function (predictors,
                 response,
                 method = "rf",
                 metric = ifelse(is.factor(response), "Accuracy", "RMSE"),
                 maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                 globalval=FALSE,
                 withinSE = FALSE,
                 minVar = 2,
                 trControl = caret::trainControl(),
                 tuneLength = 3,
                 tuneGrid = NULL,
                 seed = sample(1:1000, 1),
                 verbose=TRUE,
                 ...){
  trControl$returnResamp <- "final"
  trControl$savePredictions <- "final"
  if(inherits(response,"character")){
    response <- factor(response)
    if(metric=="RMSE"){
      metric <- "Accuracy"
      maximize <- TRUE
    }
  }
  if (trControl$method=="LOOCV"){
    if (withinSE==TRUE){
      print("warning: withinSE is set to FALSE as no SE can be calculated using method LOOCV")
      withinSE <- FALSE
    }}

  if(globalval){
    if (withinSE==TRUE){
      print("warning: withinSE is set to FALSE as no SE can be calculated using global validation")
      withinSE <- FALSE
    }}

  se <- function(x){sd(x, na.rm = TRUE)/sqrt(length(na.exclude(x)))}
  n <- length(names(predictors))

  acc <- 0
  perf_all <- data.frame(matrix(ncol=length(predictors)+3,
                                nrow=choose(n, minVar)+(n-minVar)*(n-minVar+1)/2))
  names(perf_all) <- c(paste0("var",1:length(predictors)),metric,"SE","nvar")
  if(maximize) evalfunc <- function(x){max(x,na.rm=TRUE)}
  if(!maximize) evalfunc <- function(x){min(x,na.rm=TRUE)}
  isBetter <- function (actmodelperf,bestmodelperf,
                        bestmodelperfSE=NULL,
                        maximization=FALSE,
                        withinSE=FALSE){
    if(withinSE){
      result <- ifelse (!maximization, actmodelperf < bestmodelperf-bestmodelperfSE,
                        actmodelperf > bestmodelperf+bestmodelperfSE)
    }else{
      result <- ifelse (!maximization, actmodelperf < bestmodelperf,
                        actmodelperf > bestmodelperf)
    }
    return(result)
  }
  #### chose initial best model from all combinations of two variables
  minGrid <- t(data.frame(combn(names(predictors),minVar)))
  for (i in 1:nrow(minGrid)){
    if (verbose){
      print(paste0("model using ",paste0(minGrid[i,],collapse=","), " will be trained now..." ))
    }
    set.seed(seed)
    #adaptations for pls:
    tuneGrid_orig <- tuneGrid
    tuneLength_orig <- tuneLength
    if(method=="pls"&!is.null(tuneGrid)&any(tuneGrid$ncomp>minVar)){
      tuneGrid <- data.frame(ncomp=tuneGrid[tuneGrid$ncomp<=minVar,])
      if(verbose){
        print(paste0("note: maximum ncomp is ", minVar))
      }
    }
    #adaptations for tuning of ranger:
    if(method=="ranger"&!is.null(tuneGrid)&any(tuneGrid$mtry>minVar)){
      tuneGrid$mtry <- minVar
      if(verbose){
        print("invalid value for mtry. Reset to valid range.")
      }
    }
    # adaptations for RF and minVar == 1 - tuneLength must be 1, only one mtry possible
    if(minVar==1 & method%in%c("ranger", "rf") & is.null(tuneGrid)){
      tuneLength <- minVar
    }

    #train model:
    model <- caret::train(predictors[minGrid[i,]],
                          response,
                          method=method,
                          metric=metric,
                          trControl=trControl,
                          tuneLength = tuneLength,
                          tuneGrid = tuneGrid,
                          ...)

    tuneGrid <- tuneGrid_orig
    tuneLength <- tuneLength_orig

    ### compare the model with the currently best model
    if (globalval){
      perf_stats <- global_validation(model)[names(global_validation(model))==metric]
    }else{
      perf_stats <- model$results[,names(model$results)==metric]
    }
    actmodelperf <- evalfunc(perf_stats)
    actmodelperfSE <- se(
      sapply(unique(model$resample$Resample),
             FUN=function(x){mean(model$resample[model$resample$Resample==x,
                                                 metric],na.rm=TRUE)}))
    if (i == 1){
      bestmodelperf <- actmodelperf
      bestmodelperfSE <- actmodelperfSE
      bestmodel <- model
    } else{
      if (isBetter(actmodelperf,bestmodelperf,maximization=maximize,withinSE=FALSE)){
        bestmodelperf <- actmodelperf
        bestmodelperfSE <- actmodelperfSE
        bestmodel <- model
      }
    }
    acc <- acc+1

    variablenames <- names(model$trainingData)[-length(names(model$trainingData))]
    perf_all[acc,1:length(variablenames)] <- variablenames
    perf_all[acc,(length(predictors)+1):ncol(perf_all)] <- c(actmodelperf,actmodelperfSE,length(variablenames))
    if(verbose){
      print(paste0("maximum number of models that still need to be trained: ",
                   round(choose(n, minVar)+(n-minVar)*(n-minVar+1)/2-acc,0)))
    }
  }
  #### increase the number of predictors by one (try all combinations)
  #and test if model performance increases
  selectedvars <- names(bestmodel$trainingData)[-which(
    names(bestmodel$trainingData)==".outcome")]


  if (globalval){
    selectedvars_perf <- global_validation(bestmodel)[names(global_validation(bestmodel))==metric]
  }else{
    if (maximize){
      selectedvars_perf <-max(bestmodel$results[,metric])
    }else{
      selectedvars_perf <- min(bestmodel$results[,metric])
    }
  }
  selectedvars_SE <- bestmodelperfSE
  if(verbose){
    print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
                 " with ",metric," ",round(selectedvars_perf,3)))
  }
  for (k in 1:(length(names(predictors))-minVar)){
    startvars <- names(bestmodel$trainingData)[-which(
      names(bestmodel$trainingData)==".outcome")]
    nextvars <- names(predictors)[-which(
      names(predictors)%in%startvars)]
    if (length(startvars)<(k+(minVar-1))){
      message(paste0("Note: No increase in performance found using more than ",
                     length(startvars), " variables"))
      bestmodel$selectedvars <- selectedvars
      bestmodel$selectedvars_perf <- selectedvars_perf[-length(selectedvars_perf)]
      bestmodel$selectedvars_perf_SE <- selectedvars_SE[-length(selectedvars_SE)] #!!!
      bestmodel$perf_all <- perf_all
      bestmodel$perf_all <- bestmodel$perf_all[!apply(is.na(bestmodel$perf_all), 1, all),]
      bestmodel$perf_all <- bestmodel$perf_all[colSums(!is.na(bestmodel$perf_all)) > 0]
      bestmodel$minVar <- minVar
      bestmodel$type <- "ffs"
      class(bestmodel) <- c("ffs", "train")
      return(bestmodel)
      break()
    }
    for (i in 1:length(nextvars)){
      if(verbose){
        print(paste0("model using additional variable ",nextvars[i], " will be trained now..." ))
      }
      set.seed(seed)

      #adaptation for pls:
      tuneGrid_orig <- tuneGrid
      if(method=="pls"&!is.null(tuneGrid)&any(tuneGrid$ncomp>ncol(predictors[,c(startvars,nextvars[i])]))){
        tuneGrid<- data.frame(ncomp=tuneGrid[tuneGrid$ncomp<=ncol(predictors[,c(startvars,nextvars[i])]),])
        if(verbose){
          print(paste0("note: maximum ncomp is ", ncol(predictors[,c(startvars,nextvars[i])])))
        }}
      #adaptation for ranger:
      if(method=="ranger"&!is.null(tuneGrid)&any(tuneGrid$mtry>ncol(predictors[,c(startvars,nextvars[i])]))){
        tuneGrid$mtry[tuneGrid$mtry>ncol(predictors[,c(startvars,nextvars[i])])] <- ncol(predictors[,c(startvars,nextvars[i])])
        if(verbose){
          print("invalid value for mtry. Reset to valid range.")
        }
      }

      model <- caret::train(predictors[,c(startvars,nextvars[i])],
                            response,
                            method = method,
                            metric=metric,
                            trControl = trControl,
                            tuneLength = tuneLength,
                            tuneGrid = tuneGrid,
                            ...)
      tuneGrid <- tuneGrid_orig

      if (globalval){
        perf_stats <- global_validation(model)[names(global_validation(model))==metric]
      }else{
        perf_stats <- model$results[,names(model$results)==metric]
      }
      actmodelperf <- evalfunc(perf_stats)

      actmodelperfSE <- se(
        sapply(unique(model$resample$Resample),
               FUN=function(x){mean(model$resample[model$resample$Resample==x,
                                                   metric],na.rm=TRUE)}))
      if(isBetter(actmodelperf,bestmodelperf,
                  selectedvars_SE[length(selectedvars_SE)], #SE from model with nvar-1
                  maximization=maximize,withinSE=withinSE)){
        bestmodelperf <- actmodelperf
        bestmodelperfSE <- actmodelperfSE
        bestmodel <- model
      }
      acc <- acc+1

      variablenames <- names(model$trainingData)[-length(names(model$trainingData))]
      perf_all[acc,1:length(variablenames)] <- variablenames
      perf_all[acc,(length(predictors)+1):ncol(
        perf_all)] <- c(actmodelperf,actmodelperfSE,length(variablenames))
      if(verbose){
        print(paste0("maximum number of models that still need to be trained: ",
                     round(choose(n, minVar)+(n-minVar)*(n-minVar+1)/2-acc,0)))
      }
    }
    selectedvars <- c(selectedvars,names(bestmodel$trainingData)[-which(
      names(bestmodel$trainingData)%in%c(".outcome",selectedvars))])
    selectedvars_SE <- c(selectedvars_SE,bestmodelperfSE)




    if (maximize){
      if(globalval){
        selectedvars_perf <- c(selectedvars_perf,global_validation(bestmodel)[names(global_validation(bestmodel))==metric])
      }else{
        selectedvars_perf <- c(selectedvars_perf,max(bestmodel$results[,metric]))
      }
    }
    if (!maximize){
      if(globalval){
        selectedvars_perf <- c(selectedvars_perf,global_validation(bestmodel)[names(global_validation(bestmodel))==metric])
      }else{
        selectedvars_perf <- c(selectedvars_perf,min(bestmodel$results[,metric]))
      }
    }
    if(verbose){
      print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
                   " with ",metric," ",round(selectedvars_perf[length(selectedvars_perf)],3)))
    }
  }

  # Old version that is not using global_validation:
  #    if (maximize){
  #      selectedvars_perf <- c(selectedvars_perf,max(bestmodel$results[,metric]))
  #      if(verbose){
  #        print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
  #                     " with ", metric," ",round(max(bestmodel$results[,metric]),3)))
  #      }
  #    }
  #    if (!maximize){
  #      selectedvars_perf <- c(selectedvars_perf,min(bestmodel$results[,metric]))
  #      if(verbose){
  #        print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
  #                     " with ",metric," ",round(min(bestmodel$results[,metric]),3)))
  #      }
  #    }

  bestmodel$selectedvars <- selectedvars
  bestmodel$selectedvars_perf <- selectedvars_perf
  bestmodel$selectedvars_perf_SE <- selectedvars_SE
  if(globalval){
    bestmodel$selectedvars_perf_SE <- NA
  }
  bestmodel$perf_all <- perf_all
  bestmodel$perf_all <- bestmodel$perf_all[!apply(is.na(bestmodel$perf_all), 1, all),]
  bestmodel$minVar <- minVar
  bestmodel$type <- "ffs"
  bestmodel$perf_all <- bestmodel$perf_all[colSums(!is.na(bestmodel$perf_all)) > 0]
  class(bestmodel) <- c("ffs", "train")
  return(bestmodel)
}
