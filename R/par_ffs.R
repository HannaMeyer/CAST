#' Forward Feature Selection Parallel
#' @description Parallel implementation of \code{\link{ffs}}
#'
#' @param cores Numeric. Number of CPU cores to use.
#'
#' @importFrom foreach %dopar%
#'
#'
#' @author Marvin Ludwig
#' @author Hanna Meyer
#'
#'
#' @note Only tested on Linux systems!
#'
#' @examples
#' \dontrun{
#' data(iris)
#' ffsmodel <- par_ffs(iris[,1:4],iris$Species, cores = 4)
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
#'
#' #load and prepare dataset:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' trainDat <- dat[dat$altitude==-0.3&year(dat$Date)==2012&week(dat$Date)%in%c(13:14),]
#'
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
#' ffsmodel <- par_ffs(trainDat[,predictors],trainDat$VW,method="rf",
#' tuneLength=1,trControl=ctrl, cores = 4, seed = 1)
#' ffsmodel

par_ffs <- function (predictors,
                     response,
                     method = "rf",
                     metric = ifelse(is.factor(response), "Accuracy", "RMSE"),
                     maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                     withinSE = FALSE,
                     minVar = 2,
                     trControl = caret::trainControl(),
                     tuneLength = 3,
                     tuneGrid = NULL,
                     seed = sample(1:1000, 1),
                     verbose=TRUE,
                     cores = 4,
                     ...){


  ## define functions ----
  se <- function(x){
    sd(x, na.rm = TRUE)/sqrt(length(na.exclude(x)))
  }

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



  ## initial variable selection

  initialvarsel = function(i){
    {
      if (verbose){
        print(paste0("model using ",paste0(minGrid[i,],collapse=","), " will be trained now..." ))
      }
      set.seed(seed)
      #adaptations for pls:
      tuneGrid_orig <- tuneGrid
      if(method=="pls"&!is.null(tuneGrid)&any(tuneGrid$ncomp>minVar)){
        tuneGrid <- data.frame(ncomp=tuneGrid[tuneGrid$ncomp<=minVar,])
        if(verbose){
          print(paste0("note: maximum ncomp is ", minVar))
        }
      }
      #train model:
      set.seed(seed)
      model <- caret::train(predictors[,minGrid[i,]],
                            response,
                            method=method,
                            metric=metric,
                            trControl=trControl,
                            tuneLength = tuneLength,
                            tuneGrid = tuneGrid,
                            ...)
      if(method=="pls"&!is.null(tuneGrid)&any(tuneGrid$ncomp>minVar)){
        tuneGrid <- tuneGrid_orig
      }
      actmodelperf <- evalfunc(model$results[,names(model$results)==metric])
      actmodelperfSE <- se(
        sapply(unique(model$resample$Resample),
               FUN=function(x){mean(model$resample[model$resample$Resample==x,
                                                   metric],na.rm=TRUE)}))


      variablenames <- names(model$trainingData)[-length(names(model$trainingData))]



      return(
        data.frame(t(variablenames),
                   actmodelperf,
                   actmodelperfSE,
                   length.variablenames = length(variablenames))
      )
    }
  }



  ## define static variables ---------
  n <- length(names(predictors))
  acc <- 0 # maybe unneccessary
  trControl$returnResamp <- "final"
  minGrid <- t(data.frame(combn(names(predictors),minVar)))


  ## some cleanup with if statements

  if(class(response)=="character"){
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


  if(maximize) evalfunc <- function(x){max(x,na.rm=TRUE)}
  if(!maximize) evalfunc <- function(x){min(x,na.rm=TRUE)}

  # paralell initial variables ----


  if(cores > 1){
    print(paste0("Training all initial ", minVar, " variable models on ",cores  ," cores. Verbose messages not supported."))
    doParallel::registerDoParallel(cores)

    perf_all = foreach::foreach(i = seq(nrow(minGrid)), .combine = rbind) %dopar% initialvarsel(i)
    doParallel::stopImplicitCluster()
  }else{
    perf_all = foreach::foreach(i = seq(nrow(minGrid)), .combine = rbind) %do% initialvarsel(i)
  }



  # find best initial variable combination
  bestmodelperf = ifelse(maximize, max(perf_all$actmodelperf), min(perf_all$actmodelperf))
  bestmodelperfSE = ifelse(maximize, max(perf_all$actmodelperfSE), min(perf_all$actmodelperfSE))
  bestminvars = unlist(perf_all[which(perf_all$actmodelperf == bestmodelperf),1:minVar])

  # train model with best initial variable combination again
  ## because the parallelization does not keep the model
  set.seed(seed)
  bestmodel = caret::train(predictors[,bestminvars],
                           response,
                           method=method,
                           metric=metric,
                           trControl=trControl,
                           tuneLength = tuneLength,
                           tuneGrid = tuneGrid,
                           ...)


  # end of initial combination ---------
  selectedvars <- names(bestmodel$trainingData)[-which(
    names(bestmodel$trainingData)==".outcome")]
  if (maximize){
    selectedvars_perf <- max(bestmodel$results[,metric])
  } else{
    selectedvars_perf <- min(bestmodel$results[,metric])
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

      # formating perf_all

      # reorder
      perf_all = perf_all[,c(grep("X", names(perf_all)),
                             seq(ncol(perf_all))[!seq(ncol(perf_all))%in%grep("X", names(perf_all))])]
      # rename
      names(perf_all) = gsub("X", "var", names(perf_all))
      names(perf_all)[!seq(ncol(perf_all))%in%grep("var[0-9]", names(perf_all))] = c(metric, "SE", "nvar")


      bestmodel$perf_all <- perf_all
      bestmodel$perf_all <- bestmodel$perf_all[!apply(is.na(bestmodel$perf_all), 1, all),]
      bestmodel$perf_all <- bestmodel$perf_all[colSums(!is.na(bestmodel$perf_all)) > 0]
      bestmodel$minVar <- minVar
      bestmodel$type <- "ffs"
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

      model <- caret::train(predictors[,c(startvars,nextvars[i])],
                            response,
                            method = method,
                            metric=metric,
                            trControl = trControl,
                            tuneLength = tuneLength,
                            tuneGrid = tuneGrid,
                            ...)
      tuneGrid <- tuneGrid_orig
      actmodelperf <- evalfunc(model$results[,names(model$results)==metric])
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

      perf_all = plyr::rbind.fill(
        perf_all,
        data.frame(t(variablenames), actmodelperf, actmodelperfSE, length.variablenames = length(variablenames))
      )



      if(verbose){
        print(paste0("maximum number of models that still need to be trained: ",
                     round(choose(n, minVar)+(n-minVar)*(n-minVar+1)/2-acc,0)))
      }
    }
    selectedvars <- c(selectedvars,names(bestmodel$trainingData)[-which(
      names(bestmodel$trainingData)%in%c(".outcome",selectedvars))])
    selectedvars_SE <- c(selectedvars_SE,bestmodelperfSE)
    if (maximize){
      selectedvars_perf <- c(selectedvars_perf,max(bestmodel$results[,metric]))
      if(verbose){
        print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
                     " with ", metric," ",round(max(bestmodel$results[,metric]),3)))
      }
    }
    if (!maximize){
      selectedvars_perf <- c(selectedvars_perf,min(bestmodel$results[,metric]))
      if(verbose){
        print(paste0(paste0("vars selected: ",paste(selectedvars, collapse = ',')),
                     " with ",metric," ",round(min(bestmodel$results[,metric]),3)))
      }
    }
  }

  # formating perf_all

  # reorder
  perf_all = perf_all[,c(grep("X", names(perf_all)),
                         seq(ncol(perf_all))[!seq(ncol(perf_all))%in%grep("X", names(perf_all))])]
  # rename
  names(perf_all) = gsub("X", "var", names(perf_all))
  names(perf_all)[!seq(ncol(perf_all))%in%grep("var[0-9]", names(perf_all))] = c(metric, "SE", "nvar")




  bestmodel$selectedvars <- selectedvars
  bestmodel$selectedvars_perf <- selectedvars_perf
  bestmodel$selectedvars_perf_SE <- selectedvars_SE
  bestmodel$perf_all <- perf_all
  bestmodel$perf_all <- bestmodel$perf_all[!apply(is.na(bestmodel$perf_all), 1, all),]
  bestmodel$minVar <- minVar
  bestmodel$type <- "ffs"
  bestmodel$perf_all <- bestmodel$perf_all[colSums(!is.na(bestmodel$perf_all)) > 0]


  return(bestmodel)
}


