#' Calculate Dissimilarity Index of training data
#' @description
#' This function estimates the Dissimilarity Index (DI) of
#' within the training data set used for a prediction model.
#' Predictors can be weighted based on the internal
#' variable importance of the machine learning algorithm used for model training.
#' @note
#' This function is called within \code{\link{aoa}} to estimate the DI and AOA of new data.
#' However, it may also be used on its own if only the DI of training data is of interest,
#' or to facilitate a parallelization of \code{\link{aoa}} by avoiding a repeated calculation of the DI within the training data.
#'
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds
#' @param train A data.frame containing the data used for model training. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Only required if no model is given.
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the model are used or if no model is given then of the train dataset.
#' @param CVtest list or vector. Either a list where each element contains the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point.
#' Only required if no model is given.
#' @param CVtrain list. Each element contains the data points used for training during the cross validation iteration (i.e. held back data).
#' Only required if no model is given and only required if CVtrain is not the opposite of CVtest (i.e. if a data point is not used for testing, it is used for training).
#' Relevant if some data points are excluded, e.g. when using \code{\link{nndm}}.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer.
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param LPD Logical. Indicates wheather the LPD should be calculated or not.
#'
#' @seealso \code{\link{aoa}}
#' @importFrom graphics boxplot
#' @import ggplot2
#'
#' @return A list of class \code{trainDI} containing:
#'  \item{train}{A data frame containing the training data}
#'  \item{weight}{A data frame with weights based on the variable importance.}
#'  \item{variables}{Names of the used variables}
#'  \item{catvars}{Which variables are categorial}
#'  \item{scaleparam}{Scaling parameters. Output from \code{scale}}
#'  \item{trainDist_avrg}{A data frame with the average distance of each training point to every other point}
#'  \item{trainDist_avrgmean}{The mean of trainDist_avrg. Used for normalizing the DI}
#'  \item{trainDI}{Dissimilarity Index of the training data}
#'  \item{threshold}{The DI threshold used for inside/outside AOA}
#'  \item{trainLPD}{LPD of the training data}
#'  \item{avrgLPD}{Average LPD of the training data}
#'
#'
#'
#' @export trainDI
#'
#' @author
#' Hanna Meyer, Marvin Ludwig
#'
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \doi{10.1111/2041-210X.13650}
#'
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(caret)
#' library(viridis)
#' library(ggplot2)
#'
#' # prepare sample data:
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' set.seed(100)
#' pts <- pts[1:30,]
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
#' trainDat <- extract(studyArea,pts,na.rm=FALSE)
#' trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
#'
#' # visualize data spatially:
#' plot(studyArea)
#' plot(studyArea$DEM)
#' plot(pts[,1],add=TRUE,col="black")
#'
#' # train a model:
#' set.seed(100)
#' variables <- c("DEM","NDRE.Sd","TWI")
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
#' trControl=trainControl(method="cv",number=5,savePredictions=T))
#' print(model) #note that this is a quite poor prediction model
#' prediction <- predict(studyArea,model,na.rm=TRUE)
#' plot(varImp(model,scale=FALSE))
#'
#' #...then calculate the DI of the trained model:
#' DI = trainDI(model=model)
#' plot(DI)
#'
#' # the DI can now be used to compute the AOA:
#' AOA = aoa(studyArea, model = model, trainDI = DI)
#' print(AOA)
#' plot(AOA)
#' }
#'


trainDI <- function(model = NA,
                    train = NULL,
                    variables = "all",
                    weight = NA,
                    CVtest = NULL,
                    CVtrain = NULL,
                    method="L2",
                    useWeight = TRUE,
                    LPD = FALSE){

  # get parameters if they are not provided in function call-----
  if(is.null(train)){train = aoa_get_train(model)}
  if(length(variables) == 1){
    if(variables == "all"){
      variables = aoa_get_variables(variables, model, train)
    }
  }
  if(is.na(weight)[1]){
    if(useWeight){
      weight = aoa_get_weights(model, variables = variables)
    }else{
      message("variable are not weighted. see ?aoa")
      weight <- t(data.frame(rep(1,length(variables))))
      names(weight) <- variables
    }
  }else{ #check if manually given weights are correct. otherwise ignore (set to 1):
    if(nrow(weight)!=1||ncol(weight)!=length(variables)){
      message("variable weights are not correctly specified and will be ignored. See ?aoa")
      weight <- t(data.frame(rep(1,length(variables))))
      names(weight) <- variables
    }
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }

  # get CV folds from model or from parameters
  folds <-  aoa_get_folds(model,CVtrain,CVtest)
  CVtest <- folds[[2]]
  CVtrain <- folds[[1]]

  # check for input errors -----
  if(nrow(train)<=1){stop("at least two training points need to be specified")}

  # reduce train to specified variables
  train <- train[,na.omit(match(variables, names(train)))]

  train_backup <- train

  # convert categorial variables
  catupdate <- aoa_categorial_train(train, variables, weight)

  train <- catupdate$train
  weight <- catupdate$weight

  # scale train
  train <- scale(train)

  # make sure all variables have variance
  if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
    stop("some variables in train seem to have no variance")
  }

  # save scale param for later
  scaleparam <- attributes(train)


  # multiply train data with variable weights (from variable importance)
  if(!inherits(weight, "error")&!is.null(unlist(weight))){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }


  # calculate average mean distance between training data

  trainDist_avrg <- c()
  trainDist_min <- c()

  if(method=="MD"){
    if(dim(train)[2] == 1){
      S <- matrix(stats::var(train), 1, 1)
    } else {
      S <- stats::cov(train)
    }
    S_inv <- MASS::ginv(S)
  }

  message("Computing DI of training data...")
  pb <- txtProgressBar(min = 0,
                       max = nrow(train),
                       style = 3)

  for(i in seq(nrow(train))){

    # distance to all other training data (for average)
    trainDistAll   <- .alldistfun(t(train[i,]), train,  method, S_inv=S_inv)[-1]
    trainDist_avrg <- append(trainDist_avrg, mean(trainDistAll, na.rm = TRUE))

    # calculate  distance to other training data:
    trainDist      <- matrix(.alldistfun(t(matrix(train[i,])), train, method, sorted = FALSE, S_inv))
    trainDist[i]   <- NA


    # mask of any data that are not used for training for the respective data point (using CV)
    whichfold <- NA
    if(!is.null(CVtrain)&!is.null(CVtest)){
      whichfold <-  as.numeric(which(lapply(CVtest,function(x){any(x==i)})==TRUE)) # index of the fold where i is held back
      if(length(whichfold)>1){stop("a datapoint is used for testing in more than one fold. currently this option is not implemented")}
      if(length(whichfold)!=0){ # in case that a data point is never used for testing
        trainDist[!seq(nrow(train))%in%CVtrain[[whichfold]]] <- NA # everything that is not in the training data for i is ignored
      }
      if(length(whichfold)==0){#in case that a data point is never used for testing, the distances for that point are ignored
        trainDist <- NA
      }
    }

    #######################################

    if (length(whichfold)==0){
      trainDist_min <- append(trainDist_min, NA)
    }else{
      trainDist_min <- append(trainDist_min, min(trainDist, na.rm = TRUE))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  trainDist_avrgmean <- mean(trainDist_avrg,na.rm=TRUE)



  # Dissimilarity Index of training data -----
  TrainDI <- trainDist_min/trainDist_avrgmean


  # AOA Threshold ----
  threshold_quantile <- stats::quantile(TrainDI, 0.75,na.rm=TRUE)
  threshold_iqr <- (1.5 * stats::IQR(TrainDI,na.rm=T))
  thres <- threshold_quantile + threshold_iqr
  # account for case that threshold_quantile + threshold_iqr is larger than maximum DI.
  if (thres>max(TrainDI,na.rm=T)){
    thres <- max(TrainDI,na.rm=T)
  }

  # note: previous versions of CAST derived the threshold this way:
  # thres <- grDevices::boxplot.stats(TrainDI)$stats[5]


  # calculate trainLPD and avrgLPD according to the CV folds
  if (LPD == TRUE) {
    message("Computing LPD of training data...")
    pb <- txtProgressBar(min = 0,
                         max = nrow(train),
                         style = 3)

    trainLPD <- c()
    for (j in  seq(nrow(train))) {

      # calculate  distance to other training data:
      trainDist      <- .alldistfun(t(matrix(train[j,])), train, method, sorted = FALSE, S_inv)
      DItrainDist <- trainDist/trainDist_avrgmean
      DItrainDist[j]   <- NA

      # mask of any data that are not used for training for the respective data point (using CV)
      whichfold <- NA
      if(!is.null(CVtrain)&!is.null(CVtest)){
        whichfold <- as.numeric(which(lapply(CVtest,function(x){any(x==j)})==TRUE)) # index of the fold where i is held back
        if(length(whichfold)!=0){ # in case that a data point is never used for testing
          DItrainDist[!seq(nrow(train))%in%CVtrain[[whichfold]]] <- NA # everything that is not in the training data for i is ignored
        }
        if(length(whichfold)==0){#in case that a data point is never used for testing, the distances for that point are ignored
          DItrainDist <- NA
        }
      }

      #######################################

      if (length(whichfold)==0){
        trainLPD <- append(trainLPD, NA)
      } else {
        trainLPD <- append(trainLPD, sum(DItrainDist[,1] < thres, na.rm = TRUE))
      }
      setTxtProgressBar(pb, j)
    }

    close(pb)

    # Average LPD in trainData
    avrgLPD <- round(mean(trainLPD))
  }


  # Return: trainDI Object -------

  aoa_results = list(
    train = train_backup,
    weight = weight,
    variables = variables,
    catvars = catupdate$catvars,
    scaleparam = scaleparam,
    trainDist_avrg = trainDist_avrg,
    trainDist_avrgmean = trainDist_avrgmean,
    trainDI = TrainDI,
    threshold = thres,
    method = method
  )

  if (LPD == TRUE) {
    aoa_results$trainLPD <- trainLPD
    aoa_results$avrgLPD <- avrgLPD
  }

  class(aoa_results) = "trainDI"

  return(aoa_results)
}


################################################################################
# Helper functions
################################################################################
# Encode categorial variables

aoa_categorial_train <- function(train, variables, weight){

  # get all categorial variables
  catvars <- tryCatch(names(train)[which(sapply(train[,variables], class)%in%c("factor","character"))],
                      error=function(e) e)

  if (!inherits(catvars,"error")&length(catvars)>0){
    message("warning: predictors contain categorical variables. The integration is currently still under development. Please check results carefully!")

    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      train[,catvar]<-droplevels(train[,catvar])

      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = train),
                           train)
      train <- data.frame(train,dvi_train)

      if(!inherits(weight, "error")){
        addweights <- data.frame(t(rep(weight[,which(names(weight)==catvar)],
                                       ncol(dvi_train))))
        names(addweights)<- colnames(dvi_train)
        weight <- data.frame(weight,addweights)
      }
    }
    if(!inherits(weight, "error")){
      weight <- weight[,-which(names(weight)%in%catvars)]
    }
    train <- train[,-which(names(train)%in%catvars)]
  }
  return(list(train = train, weight = weight, catvars = catvars))


}



# Get weights from train object


aoa_get_weights = function(model, variables){

  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error") & length(variables)>1){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    # set all weights to 1
    weight <- as.data.frame(t(rep(1, length(variables))))
    names(weight) = variables
    message("note: variables were not weighted either because no weights or model were given,
    no variable importance could be retrieved from the given model, or the model has a single feature.
    Check caret::varImp(model)")
  }

  #set negative weights to 0
  if(!inherits(weight, "error")){
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  return(weight)

}


# Get trainingdata from train object

aoa_get_train <- function(model){

  train <- as.data.frame(model$trainingData)
  return(train)


}


# Get folds from train object


aoa_get_folds <- function(model, CVtrain, CVtest){
  ### if folds are to be extracted from the model:
  if (!is.na(model)[1]){
    if(tolower(model$control$method)!="cv"){
      message("note: Either no model was given or no CV was used for model training. The DI threshold is therefore based on all training data")
    }else{
      CVtest <- model$control$indexOut
      CVtrain <- model$control$index
    }
  }
  ### if folds are specified manually:
  if(is.na(model)[1]){

    if(!is.null(CVtest)&!is.list(CVtest)){ # restructure input if CVtest only contains the fold ID
      tmp <- list()
      for (i in unique(CVtest)){
        tmp[[i]] <- which(CVtest==i)
      }
      CVtest <- tmp
    }

    if(is.null(CVtest)&is.null(CVtrain)){
      message("note: No model and no CV folds were given. The DI threshold is therefore based on all training data")
    }else{
      if(is.null(CVtest)){ # if CVtest is not given, then use the opposite of CVtrain
        CVtest <- lapply(CVtrain,function(x){which(!sort(unique(unlist(CVtrain)))%in%x)})
      }else{
        if(is.null(CVtrain)){ # if CVtrain is not given, then use the opposite of CVtest
          CVtrain <- lapply(CVtest,function(x){which(!sort(unique(unlist(CVtest)))%in%x)})
        }
      }
    }
  }
  return(list(CVtrain,CVtest))
}






# Get variables from train object

aoa_get_variables <- function(variables, model, train){

  if(length(variables) == 1){
    if(variables == "all"){
      if(!is.na(model)[1]){
        variables <- names(model$trainingData)[-which(names(model$trainingData)==".outcome")]
      }else{
        variables <- names(train)
      }
    }
  }
  return(variables)


}



.mindistfun <- function(point, reference, method, S_inv=NULL){

  if (method == "L2"){ # Euclidean Distance
    return(c(FNN::knnx.dist(reference, point, k = 1)))
  } else if (method == "MD"){ # Mahalanobis Distance
    return(sapply(1:dim(point)[1],
                  function(y) min(sapply(1:dim(reference)[1],
                                         function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) )))))
  }
}

.alldistfun <- function(point, reference, method, sorted = TRUE,S_inv=NULL){

  if (method == "L2"){ # Euclidean Distance
    if(sorted){
      return(FNN::knnx.dist(reference, point, k = dim(reference)[1]))
    } else {
      return(FNN::knnx.dist(point,reference,k=1))
    }
  } else if (method == "MD"){ # Mahalanobis Distance
    if(sorted){
      return(t(sapply(1:dim(point)[1],
                      function(y) sort(sapply(1:dim(reference)[1],
                                              function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) ))))))
    } else {
      return(t(sapply(1:dim(point)[1],
                      function(y) sapply(1:dim(reference)[1],
                                         function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) )))))
    }
  }
}

