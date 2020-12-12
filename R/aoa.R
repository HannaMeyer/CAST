#' Area of Applicability
#' @description
#' This function estimates the Dissimilarity Index (DI) and the derived
#' Area of Applicability (AOA) of spatial prediction models by
#' considering the distance of new data (i.e. a Raster Stack of spatial predictors
#' used in the models) in the predictor variable space to the data used for model
#' training. Predictors can be weighted in the ideal case based on the internal
#' variable importance of the machine learning algorithm used for model training.
#' @param newdata A RasterStack, RasterBrick or data.frame containing the data
#' the model was meant to make predictions for.
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds
#' @param cl A cluster object e.g. created with doParallel. Should only be used if newdata is large.
#' @param train A data.frame containing the data used for model training. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Only required if no model is given.
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the model are used or if no model is given then of the train dataset.
#' @param folds Numeric or character. Folds for cross validation. E.g. Spatial cluster affiliation for each data point.
#' Should be used if replicates are present. Only required if no model is given.
#' @param returnTrainDI A logical: should the DI value of the cross-validated training data be returned as a attribute?
#' @details The Dissimilarity Index (DI) and the corresponding Area of Applicability (AOA) are calculated.
#' If variables are factors, dummy variables are created prior to weighting and distance calculation.
#'
#' Interpretation of results: If a location is very similar to the properties
#' of the training data it will have a low distance in the predictor variable space
#' (DI towards 0) while locations that are very different in their properties
#' will have a high DI.
#' To get the AOA, a threshold to the DI is applied based on the DI in the training data.
#' To calculate the DI in the training data, the minimum distance to an other training point
#' (if applicable: not located in the same CV fold) is considered.
#' The DI threshold is 1,5×IQR of the DI of the training data.
#' See Meyer and Pebesma (2020) for the full documentation of the methodology.
#' @return A RasterStack or data.frame with the DI and AOA. AOA has values 0 (outside AOA) and 1 (inside AOA).
#' @author
#' Hanna Meyer
#' @references Meyer, H., Pebesma, E. (2020): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \url{https://arxiv.org/abs/2005.07939}
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#' library(viridis)
#' library(latticeExtra)
#'
#' # prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' set.seed(100)
#' pts <- pts[1:30,]
#' studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))[[1:8]]
#' trainDat <- extract(studyArea,pts,df=TRUE)
#' trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
#'
#' # visualize data spatially:
#' spplot(scale(studyArea))
#' plot(studyArea$DEM)
#' plot(pts[,1],add=TRUE,col="black")
#'
#' # train a model:
#' set.seed(100)
#' variables <- c("DEM","NDRE.Sd","TWI")
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW,method="rf",importance=TRUE,tuneLength=1,trControl=trainControl(method="cv",number=5,savePredictions=T))
#' print(model) #note that this is a quite poor prediction model
#' prediction <- predict(studyArea,model)
#' plot(varImp(model,scale=FALSE))
#'
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea,model)
#' spplot(AOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' #plot predictions for the AOA only:
#' spplot(prediction, col.regions=viridis(100),main="prediction for the AOA")+
#' spplot(AOA$AOA,col.regions=c("grey","transparent"))
#'
#' ####
#' # Calculating the AOA might be time consuming. Consider running it in parallel:
#' ####
#' library(doParallel)
#' library(parallel)
#' cl <- makeCluster(4)
#' registerDoParallel(cl)
#' AOA <- aoa(studyArea,model,cl=cl)
#'
#' ####
#' #The AOA can also be calculated without a trained model.
#' #All variables are weighted equally in this case:
#' ####
#' AOA <- aoa(studyArea,train=trainDat,variables=variables)
#' spplot(AOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' spplot(AOA$AOA,main="Area of Applicability")
#' }
#' @export aoa
#' @aliases aoa

aoa <- function(newdata,
                model=NA,
                cl=NULL,
                train=NULL,
                weight=NA,
                variables="all",
                folds=NULL,
                returnTrainDI=TRUE) {

  ### if not specified take all variables from train dataset as default:
  if(is.null(train)){train <- model$trainingData}
  if(nrow(train)<=1){stop("at least two training points need to be specified")}
  if(length(variables)==1&&variables=="all"){
    if(!is.na(model)[1]){
      variables <- names(model$trainingData)[-which(names(model$trainingData)==".outcome")]
    }else{
      variables <- names(train)
    }
  }
  as_stars <- FALSE
  if (inherits(newdata, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    newdata = as(newdata, "Raster")
    as_stars <- TRUE
  }
  #### Prepare output as either as RasterLayer or vector:
  out <- NA
  if (inherits(newdata, "Raster"))
    out <- newdata[[1]]

  #### Extract weights from trained model:
  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error")){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    message("note: variables were not weighted either because no weights or model were given,
    or no variable importance could be retrieved from the given model.
    Check caret::varImp(model)")
  }

  #### order data:
  if (inherits(newdata, "Raster")){
    if (any(is.factor(newdata))){
      newdata[[which(is.factor(newdata))]] <- raster::deratify(newdata[[which(is.factor(newdata))]],complete = TRUE)
    }
    newdata <- raster::as.data.frame(newdata)
  }
  newdata <- newdata[,na.omit(match(variables, names(newdata)))]
  train <- train[,na.omit(match(variables, names(train)))]
  if(!inherits(weight, "error")){
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  ##############################################################################
  ## Handling of categorical predictors:
  catvars <- tryCatch(names(train)[which(sapply(train[,variables], class)%in%c("factor","character"))], error=function(e) e)
  if (!inherits(catvars,"error")&length(catvars)>0){
    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      train[,catvar]<-droplevels(train[,catvar])
      newdata[,catvar] <- factor(newdata[,catvar])
      newdata[!newdata[,catvar]%in%unique(train[,catvar]),catvar] <- NA
      newdata[,catvar] <- droplevels(newdata[,catvar])
      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = train),train)
      dvi_newdata <- predict(caret::dummyVars(paste0("~",catvar), data=train),newdata)
      dvi_newdata[is.na(newdata[,catvar]),] <- 0
      train <- data.frame(train,dvi_train)
      newdata <- data.frame(newdata,dvi_newdata)
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
    newdata <- newdata[,-which(names(newdata)%in%catvars)]
    train <- train[,-which(names(train)%in%catvars)]
  }
  ##############################################################################
  #### Scale data and weight data if applicable:
  train <- scale(train)
  if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
    stop("some variables in train seem to have no variance")
  }
  scaleparam <- attributes(train)
  if(!inherits(weight, "error")){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }
  newdata <- scale(newdata,center=scaleparam$`scaled:center`,#scaleparam$`scaled:center`
                   scale=scaleparam$`scaled:scale`)

  if(!inherits(weight, "error")){
    newdata <- sapply(1:ncol(newdata),function(x){newdata[,x]*unlist(weight[x])})
  }

  #### For each pixel caclculate distance to each training point and search for
  #### min distance:

  distfun <- function(x){
    if(any(is.na(x))){
      return(NA)
    }else{
      tmp <- FNN::knnx.dist(t(matrix(x)),train,k=1)
      return(min(tmp))
    }
  }
  if (!is.null(cl)){
    mindist <- parApply(cl=cl,X=newdata,MARGIN=1,FUN=distfun)
  }else{
    mindist <- apply(newdata,1,FUN=distfun)
  }

  if(!is.null(model)&is.null(folds)){
    CVfolds <- tryCatch(reshape::melt(model$control$indexOut),
                        error=function(e) e)
    if(inherits(CVfolds, "error")){
      message("note: Either no model was given or no CV was used for model training. The DI threshold is therefore based on all training data")
    }else{
      CVfolds <- CVfolds[order(CVfolds$value),]
    }
  }

  # calculate average mean distance between training data

  trainDist_avrg <-  apply(train,1,FUN=function(x){
    mean(FNN::knnx.dist(train, t(data.frame(x)),
                   k=nrow(train))[-1],na.rm=TRUE)

  })
  trainDist_avrgmean <- mean(trainDist_avrg,na.rm=TRUE)


 # trainDist_mean <- c()
  trainDist_min <- c()
  for (i in 1:nrow(train)){
    #calculate distance to other training data:
    trainDist <- FNN::knnx.dist(t(matrix(train[i,])),train,k=1)
    trainDist[i] <- NA

    # If data are highly clustered (repliates) make sure that distance to data from same
    # cluster are excluded. Only required if no model with CV folds is given:
    if (!is.null(folds)){
      trainDist[folds==folds[i]] <- NA
    }

    # if folds are not manually assigned, CV folds from the model will be used
    # to derive the threshold on the DI:
    if(!is.na(model[1])){
      trainDist[CVfolds$L1==CVfolds$L1[i]] <- NA
    }
    # get minimum and mean distance to other training locations:
   # trainDist_mean <- c(trainDist_mean,mean(trainDist,na.rm=T))
    trainDist_min <- c(trainDist_min,min(trainDist,na.rm=T))
  }

  #scale the distance to nearest training point by average distance of the training data
  #trainDist_avrgmean <- mean(trainDist_mean)
  DI_out <- mindist/trainDist_avrgmean

  # define threshold for AOA:
  TrainDI <- trainDist_min/trainDist_avrgmean
  AOA_train_stats <- quantile(TrainDI,
                              probs = c(0.25,0.5,0.75,0.9,0.95,0.99,1),na.rm = TRUE)
#  if(is.null(thres)){
    thres <- boxplot.stats(TrainDI)$stats[5]
    lower_thres <- boxplot.stats(TrainDI)$stats[1]
#  }else{
#    thres <- quantile(TrainDI,probs = thres,na.rm=TRUE)
#  }
  #### Create Mask for AOA and return statistics
  if (inherits(out, "RasterLayer")){
    raster::values(out) <- DI_out

    AOA <- out
    raster::values(AOA) <- 1
    AOA[out>thres] <- 0
    AOA <- raster::mask(AOA,out)

    CVA <- out
    raster::values(CVA) <- 1
    CVA[out>thres|out<lower_thres] <- 0
    CVA <- raster::mask(CVA,out)


    out <- raster::stack(out,AOA,CVA)
    if (as_stars)
      out <- split(stars::st_as_stars(out), "band")
  }else{
    out <- DI_out
    AOA <- rep(1,length(out))
    AOA[out>thres] <- 0

    CVA <- rep(1,length(out))
    CVA[out>thres|out<lower_thres] <- 0

    out <- list(out,AOA,CVA)
  }
  names(out) <- c("DI","AOA","CVA")
  attributes(out)$aoa_stats <- list("Mean_train" = trainDist_avrgmean,
                                    "threshold_stats" = AOA_train_stats,
                                    "threshold" = thres,
                                    "lower_threshold"=lower_thres)
  if(returnTrainDI){
    attributes(out)$TrainDI <- TrainDI
  }
  return(out)
}
