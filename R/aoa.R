#' Area of Applicability
#' @description
#' This function estimates the Dissimilarity Index (DI) and the derived
#' Area of Applicability (AOA) of spatial prediction models by
#' considering the distance of new data (i.e. a Raster Stack of spatial predictors
#' used in the models) in the predictor variable space to the data used for model
#' training. Predictors can be weighted based on the internal
#' variable importance of the machine learning algorithm used for model training.
#' The AOA is derived by applying a threshold on the DI which is the (outlier-removed)
#' maximum DI of the cross-validated training data.
#' @param newdata A RasterStack, RasterBrick, stars object, SpatRaster or data.frame containing the data
#' the model was meant to make predictions for.
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds
#' @param trainDI A trainDI object. Optional if \code{\link{trainDI}} was calculated beforehand.
#' @param cl A cluster object e.g. created with doParallel. Should only be used if newdata is large.
#' @param train A data.frame containing the data used for model training. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Only required if no model is given.
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the model are used or if no model is given then of the train dataset.
#' @param folds Numeric or character. Folds for cross validation. E.g. Spatial cluster affiliation for each data point.
#' Should be used if replicates are present. Only required if no model is given.
#' @details The Dissimilarity Index (DI) and the corresponding Area of Applicability (AOA) are calculated.
#' If variables are factors, dummy variables are created prior to weighting and distance calculation.
#'
#' Interpretation of results: If a location is very similar to the properties
#' of the training data it will have a low distance in the predictor variable space
#' (DI towards 0) while locations that are very different in their properties
#' will have a high DI.
#' See Meyer and Pebesma (2020) for the full documentation of the methodology.
#' @note If classification models are used, currently the variable importance can only
#' be automatically retrieved if models were traiend via train(predictors,response) and not via the formula-interface.
#' Will be fixed.
#' @return A list of class \code{aoa} containing:
#'  \item{parameters}{trainDI parameters. see \code{\link{trainDI}}}
#'  \item{DI}{raster or data frame. Dissimilarity index of newdata}
#'  \item{AOA}{raster or data frame. Area of Applicability of newdata.
#'   AOA has values 0 (outside AOA) and 1 (inside AOA)}
#'
#' @author
#' Hanna Meyer
#' @references Meyer, H., Pebesma, E. (2020): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \url{https://arxiv.org/abs/2005.07939}
#' @seealso \code{\link{calibrate_aoa}}
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
#' trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
#' trControl=trainControl(method="cv",number=5,savePredictions=T))
#' print(model) #note that this is a quite poor prediction model
#' prediction <- predict(studyArea,model)
#' plot(varImp(model,scale=FALSE))
#'
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea,model)
#' spplot(AOA$predictionAOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' #plot predictions for the AOA only:
#' spplot(prediction, col.regions=viridis(100),main="prediction for the AOA")+
#' spplot(AOA$predictionAOA$AOA,col.regions=c("grey","transparent"))
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
#' spplot(AOA$predictionAOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' spplot(AOA$predictionAOA$AOA,main="Area of Applicability")
#' }
#' @export aoa
#' @aliases aoa

aoa <- function(newdata,
                model=NA,
                trainDI = NA,
                cl=NULL,
                train=NULL,
                weight=NA,
                variables="all",
                folds=NULL) {

  # if not provided, compute train DI
  if(!inherits(trainDI, "trainDI")){
    message("No trainDI provided.")
    message("Computing DI of training data...")
    trainDI = trainDI(model, train, variables, weight, folds)
  }

  message("Computing DI of newdata...")


  # check if variables are in newdata
  if(any(trainDI$variables %in% names(newdata)==FALSE)){
    stop("names of newdata don't match names of train data in the model")
  }


  # handling of different raster formats
  as_stars <- FALSE
  as_terra <- FALSE
  if (inherits(newdata, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    newdata = methods::as(newdata, "Raster")
    as_stars <- TRUE
  }

  if (inherits(newdata, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("package terra required: install that first")
    newdata = methods::as(newdata, "Raster")
    as_terra <- TRUE
  }



  # Prepare output as either as RasterLayer or vector:
  out <- NA
  if (inherits(newdata, "Raster")){
    out <- newdata[[1]]
  }



  #### order data:
  if (inherits(newdata, "Raster")){
    if (any(is.factor(newdata))){
      newdata[[which(is.factor(newdata))]] <- raster::deratify(newdata[[which(is.factor(newdata))]],complete = TRUE)
    }
    newdata <- raster::as.data.frame(newdata)
  }
  newdata <- newdata[,na.omit(match(trainDI$variables, names(newdata)))]



  ## Handling of categorical predictors:
  catvars <- trainDI$catvars
  if (!inherits(catvars,"error")&length(catvars)>0){
    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      trainDI$train[,catvar]<-droplevels(trainDI$train[,catvar])
      newdata[,catvar] <- factor(newdata[,catvar])
      newdata[!newdata[,catvar]%in%unique(trainDI$train[,catvar]),catvar] <- NA
      newdata[,catvar] <- droplevels(newdata[,catvar])
      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = trainDI$train),trainDI$train)
      dvi_newdata <- predict(caret::dummyVars(paste0("~",catvar), data=trainDI$train),newdata)
      dvi_newdata[is.na(newdata[,catvar]),] <- 0
      trainDI$train <- data.frame(trainDI$train,dvi_train)
      newdata <- data.frame(newdata,dvi_newdata)

    }
    newdata <- newdata[,-which(names(newdata)%in%catvars)]
    trainDI$train <- trainDI$train[,-which(names(trainDI$train)%in%catvars)]
  }


  newdata <- scale(newdata,center=trainDI$scaleparam$`scaled:center`,
                   scale=trainDI$scaleparam$`scaled:scale`)

  if(!inherits(trainDI$weight, "error")){
    newdata <- sapply(1:ncol(newdata),function(x){newdata[,x]*unlist(trainDI$weight[x])})
  }



  # Distance Calculation ---------


  distfun <- function(x){
    if(any(is.na(x))){
      return(NA)
    }else{

      # rescale and reweight train data
      train_scaled = scale(trainDI$train,
                           center = trainDI$scaleparam$`scaled:center`,
                           scale = trainDI$scaleparam$`scaled:scale`)

      train_scaled <- sapply(1:ncol(train_scaled),function(x){train_scaled[,x]*unlist(trainDI$weight[x])})

      tmp <- FNN::knnx.dist(t(matrix(x)), train_scaled, k=1)
      return(min(tmp))
    }
  }
  if (!is.null(cl)){
    mindist <- parallel::parApply(cl=cl,X=newdata,MARGIN=1,FUN=distfun)
  }else{
    mindist <- apply(newdata,1,FUN=distfun)
  }



  DI_out <- mindist/trainDI$trainDist_avrgmean


  message("Computing AOA...")

  #### Create Mask for AOA and return statistics
  if (inherits(out, "RasterLayer")){
    raster::values(out) <- DI_out

    AOA <- out
    raster::values(AOA) <- 1
    AOA[out>trainDI$thres] <- 0
    AOA <- raster::mask(AOA,out)

    #CVA <- out
    #raster::values(CVA) <- 1
    #CVA[out>thres|out<lower_thres] <- 0
    #CVA <- raster::mask(CVA,out)
    #out <- raster::stack(out,AOA,CVA)





    # handling of different raster formats. Eventually outsource to full method.

    if (as_stars){
      out <- stars::st_as_stars(out)
      AOA <- stars::st_as_stars(AOA)
    }

    if(as_terra){
      out = methods::as(out, "SpatRaster")
      AOA = methods::as(AOA, "SpatRaster")
    }




  }else{
    out <- data.frame(DI = DI_out)
    AOA <- rep(1,nrow(out))
    AOA[out>trainDI$thres] <- 0
    AOA <- as.data.frame(AOA)

    #CVA <- rep(1,length(out))
    #CVA[out>trainDI$thres|out<trainDI$lower_thres] <- 0

    #out <- list(out,AOA,CVA)

  }


  # eventually remove this
  attributes(AOA)$aoa_stats <- list("Mean_train" = trainDI$trainDist_avrgmean,
                                    "threshold_stats" = trainDI$AOA_train_stats,
                                    "threshold" = trainDI$thres,
                                    "lower_threshold" = trainDI$lower_thres)
  attributes(AOA)$TrainDI = trainDI$trainDI




  outout = list(parameters = trainDI,
                DI = out,
                AOA = AOA)

  class(outout) = "aoa"
  return(outout)


}
