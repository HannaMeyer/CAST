#' Estimate uncertainty of spatial prediction models
#'
#' @description
#' this function estimates uncertainty of spatial prediction models by
#' considering the distance of new data (i.e. a Raster Stack of spatial predictors
#' used in the models) in the predictor variable space to the data used for model
#' training. Predictors can be weighted in the ideal case based on the internal
#' variable importance of the machine learning algorithm used for model training.
#'
#' @param train a data.frame containing the data used for model training
#' @param predictors A RasterStack, RasterBrick or data.frame containing the data
#' the model was meant to make predictions for.
#' @param weights A data.frame containing weights for each variable
#' @param model A caret model used to extract weights from (based on variable importance)
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the train dataset are used.
#' @details Results range from 0 to 1 with 0 has only a short distance to the nearest
#' training point and predictions can be considered very certain. 1 has the largest
#' distance to the training data, hence predictions are more uncertain.
#' Or explained in an other way: If a location is very similar to the properties
#' of the training data it will have a low uncertainty while locations that are
#' very different in its properties will have a high uncertainty.
#' @return A RasterLayer or data.frame containing the scaled distance to the
#' nearest training data point in the predictor space.
#' @author
#' Hanna Meyer
#'
#' @examples
#' library(sf)
#' library(raster)
#' library(caret)
#' library(lubridate)
#'
#' # prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#' variables <- c("DEM","Easting","Northing")
#' trainDat <- aggregate(dat[,c("VW",variables)],by=list(as.character(dat$SOURCEID)),mean)
#'
#' # visualize data spatially:
#' plot(studyArea[[1]])
#' pts <- st_as_sf(trainDat,coords=c("Easting","Northing"))
#' plot(pts["Group.1"],add=TRUE,col="black")
#'
#' # first calculate uncertainty based on a set of variables with equal weights:
#' plot(uncertainty(trainDat,studyArea,variables=variables))
#' plot(pts["Group.1"],add=TRUE,col="black") #add training data to plot
#'
#' # or weight variables based on variable improtance from a trained model:
#' set.seed(100)
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW,method="rf",importance=TRUE,tuneLength=1)
#' plot(varImp(model))
#' # note that coordinates are the major predictors here, so uncertainty becomes higher
#' # when moving away from the training data:
#' plot(uncertainty(trainDat,studyArea,model=model,variables=variables))
#' plot(pts["Group.1"],add=TRUE,col="black") #add training data to plot
#'
#'
#' @export uncertainty
#' @aliases uncertainty

uncertainty <- function (train, predictors, weights=NA, model=NA,
                         variables="all"){
  ### if not specified take all variables from train dataset as default:
  if(length(variables)==1&&variables=="all"){
    variables=names(train)
  }
  #### Prepare output as either as RasterLayer or vector:
  out <- NA
  if (class(predictors)=="RasterStack"|class(predictors)=="RasterBrick"){
    out <- predictors[[1]]
  }
  #### Extract weights from trained model:
  #!!! achtung hier noch nach overall schauen!!!! statt [,1] für classifications
  weights <- tryCatch(as.data.frame(t(caret::varImp(model)$importance[,1])),
                      error=function(e) e)
  if(!inherits(weights, "error")){
    names(weights)<- rownames(caret::varImp(model)$importance)
  }
  #### order data:
  predictors <- predictors[[na.omit(match(variables, names(predictors)))]]
  train <- train[,na.omit(match(variables, names(train)))]
  if(!inherits(weights, "error")){
    weights <- weights[,na.omit(match(variables, names(weights)))]
  }
  #### Scale data and weight predictors if applicable:
  train <- scale(train)
  scaleparam <- attributes(train)
  if(!inherits(weights, "error")){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weights[x])})
  }
  if (class(predictors)=="RasterStack"|class(predictors)=="RasterBrick"){
    predictors <- raster::as.data.frame(predictors)
  }
  predictors <- scale(predictors,center=scaleparam$`scaled:center`,
                      scale=scaleparam$`scaled:scale`)
  if(!inherits(weights, "error")){
    predictors <- sapply(1:ncol(predictors),function(x){predictors[,x]*unlist(weights[x])})
  }
  #### For each pixel caclculate distance to each training point and search for
  #### min distance:
  tmp <- NA
  for (i in 1:nrow(train)){
    mindist <- apply(predictors,1,function(x){dist(rbind(x,train[i,]))})
    mindist <- pmin(mindist,tmp,na.rm=T)
    tmp <- mindist
  }
  #### return scaled distances as RasterLayer or vector:
  if (class(out)=="RasterLayer"){
    raster::values(out) <- scales::rescale(mindist, to = c(0, 1))
  } else{
    out <- scales::rescale(mindist, to = c(0, 1))
  }
  return(out)
}

