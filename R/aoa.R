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
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds.
#' See examples for the case that no model is available or for models trained via e.g. mlr3.
#' @param trainDI A trainDI object. Optional if \code{\link{trainDI}} was calculated beforehand.
#' @param train A data.frame containing the data used for model training. Optional. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Optional. Only required if no model is given.
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
#' @details The Dissimilarity Index (DI) and the corresponding Area of Applicability (AOA) are calculated.
#' If variables are factors, dummy variables are created prior to weighting and distance calculation.
#'
#' Interpretation of results: If a location is very similar to the properties
#' of the training data it will have a low distance in the predictor variable space
#' (DI towards 0) while locations that are very different in their properties
#' will have a high DI.
#' See Meyer and Pebesma (2021) for the full documentation of the methodology.
#' @note If classification models are used, currently the variable importance can only
#' be automatically retrieved if models were trained via train(predictors,response) and not via the formula-interface.
#' Will be fixed.
#' @return An object of class \code{aoa} containing:
#'  \item{parameters}{object of class trainDI. see \code{\link{trainDI}}}
#'  \item{DI}{raster or data frame. Dissimilarity index of newdata}
#'  \item{AOA}{raster or data frame. Area of Applicability of newdata.
#'   AOA has values 0 (outside AOA) and 1 (inside AOA)}
#'
#' @author
#' Hanna Meyer
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' Methods in Ecology and Evolution 12: 1620-1633. \doi{10.1111/2041-210X.13650}
#' @seealso \code{\link{calibrate_aoa}}, \code{\link{trainDI}}
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
#' plot(AOA)
#' spplot(AOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' #plot predictions for the AOA only:
#' spplot(prediction, col.regions=viridis(100),main="prediction for the AOA")+
#' spplot(AOA$AOA,col.regions=c("grey","transparent"))
#'
#' ####
#' #The AOA can also be calculated without a trained model.
#' #All variables are weighted equally in this case:
#' ####
#' AOA <- aoa(studyArea,train=trainDat,variables=variables)
#' spplot(AOA$DI, col.regions=viridis(100),main="Dissimilarity Index")
#' spplot(AOA$AOA,main="Area of Applicability")
#'
#'
#' ####
#' # The AOA can also be used for models trained via mlr3 (parameters have to be assigned manually):
#' ####
#'
#' library(mlr3)
#' library(mlr3learners)
#' library(mlr3spatial)
#' library(mlr3spatiotempcv)
#' library(mlr3extralearners)
#'
#' # initiate and train model:
#' train_df <- trainDat[, c("DEM","NDRE.Sd","TWI", "VW")]
#' backend <- as_data_backend(train_df)
#' task <- as_task_regr(backend, target = "VW")
#' lrn <- lrn("regr.randomForest", importance = "mse")
#' lrn$train(task)
#'
#' # cross-validation folds
#' rsmp_cv <- rsmp("cv", folds = 5L)$instantiate(task)
#'
#' ## predict:
#' prediction <- predict(studyArea,lrn$model)
#'
#' ### Estimate AOA
#' AOA <- aoa(studyArea,
#'            train = as.data.frame(task$data()),
#'            variables = task$feature_names,
#'            weight = data.frame(t(lrn$importance())),
#'            CVtest = rsmp_cv$instance[order(row_id)]$fold)
#'
#' }
#' @export aoa
#' @aliases aoa

aoa <- function(newdata,
                model=NA,
                trainDI = NA,
                train=NULL,
                weight=NA,
                variables="all",
                CVtest=NULL,
                CVtrain=NULL,
                method="L2",
                useWeight=TRUE) {

  # handling of different raster formats
  as_stars <- FALSE
  as_terra <- FALSE
  leading_digit <- all(grepl("[[:digit:]]",names(newdata)))

  if (inherits(newdata, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    newdata = methods::as(newdata, "Raster")
    as_stars <- TRUE
  }

  if (inherits(newdata, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("package terra required: install that first")

    newdata <- methods::as(newdata, "Raster")
    as_terra <- TRUE
  }




  # if not provided, compute train DI
  if(!inherits(trainDI, "trainDI")){
    message("No trainDI provided. Computing DI of training data...")
    trainDI <- trainDI(model, train, variables, weight, CVtest, CVtrain,method, useWeight)
  }

  message("Computing DI of newdata...")


  # check if variables are in newdata
  if(any(trainDI$variables %in% names(newdata)==FALSE)){
    if(!leading_digit)
      stop("names of newdata start with leading digits, automatically added 'X' results in mismatching names of train data in the model")
    stop("names of newdata don't match names of train data in the model")
  }


  # Prepare output as either as RasterLayer or vector:
  out <- NA
  if (inherits(newdata, "Raster")){
    out <- newdata[[1]]
    names(out) <- "DI"
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

  # scale and weight new data
  newdata <- scale(newdata,center=trainDI$scaleparam$`scaled:center`,
                   scale=trainDI$scaleparam$`scaled:scale`)

  if(!inherits(trainDI$weight, "error")){
    newdata <- sapply(1:ncol(newdata),function(x){
      newdata[,x]*unlist(trainDI$weight[x])
      })
  }


  # rescale and reweight train data
  train_scaled <- scale(trainDI$train,
                       center = trainDI$scaleparam$`scaled:center`,
                       scale = trainDI$scaleparam$`scaled:scale`)

  train_scaled <- sapply(1:ncol(train_scaled),function(x){train_scaled[,x]*unlist(trainDI$weight[x])})


  # Distance Calculation ---------

  mindist <- rep(NA, nrow(newdata))
  okrows <- which(apply(newdata, 1, function(x) all(!is.na(x))))
  newdataCC <- newdata[okrows,]


  if(method=="MD"){
    if(dim(train_scaled)[2] == 1){
      S <- matrix(stats::var(train_scaled), 1, 1)
      newdataCC <- as.matrix(newdataCC,ncol=1)
    } else {
      S <- stats::cov(train_scaled)
    }
    S_inv <- MASS::ginv(S)
  }

  mindist[okrows] <- .mindistfun(newdataCC, train_scaled, method, S_inv)


  DI_out <- mindist/trainDI$trainDist_avrgmean


  message("Computing AOA...")

  #### Create Mask for AOA and return statistics
  if (inherits(out, "RasterLayer")){
    raster::values(out) <- DI_out
    AOA <- out
    raster::values(AOA) <- 1
    AOA[out>trainDI$thres] <- 0
    AOA <- raster::mask(AOA,out)
    names(AOA) = "AOA"


    # handling of different raster formats.
    if (as_stars){
      out <- stars::st_as_stars(out)
      AOA <- stars::st_as_stars(AOA)
    }
    if(as_terra){
      out <- methods::as(out, "SpatRaster")
      AOA <- methods::as(AOA, "SpatRaster")
    }

  }else{
    out <- DI_out
    AOA <- rep(1,length(out))
    AOA[out>trainDI$thres] <- 0
  }


  # used in old versions of the AOA. eventually remove the attributes
  attributes(AOA)$aoa_stats <- list("Mean_train" = trainDI$trainDist_avrgmean,
                                    "threshold" = trainDI$thres)
  attributes(AOA)$TrainDI <- trainDI$trainDI

  result <- list(parameters = trainDI,
                DI = out,
                AOA = AOA)

  class(result) <- "aoa"
  return(result)

}


