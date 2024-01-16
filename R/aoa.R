#' Area of Applicability
#' @description
#' This function estimates the Dissimilarity Index (DI) and the derived
#' Area of Applicability (AOA) of spatial prediction models by
#' considering the distance of new data (i.e. a SpatRaster of spatial predictors
#' used in the models) in the predictor variable space to the data used for model
#' training. Predictors can be weighted based on the internal
#' variable importance of the machine learning algorithm used for model training.
#' The AOA is derived by applying a threshold on the DI which is the (outlier-removed)
#' maximum DI of the cross-validated training data.
#' @param newdata A SpatRaster, stars object or data.frame containing the data
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
#' @param LPD Logical. Indicates wheather the LPD should be calculated or not.
#' @param maxLPD numeric or integer. Only if \code{LPD = TRUE}. Number of nearest neighbors to be considered for the calculation of the LPD. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples. CAUTION! If not all training samples are considered, a fitted relationship between LPD and error metric will not make sense (@seealso \code{\link{LPDtoErrormetri}})
#' @details The Dissimilarity Index (DI), the Local Data Point Density (LPD) and the corresponding Area of Applicability (AOA) are calculated.
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
#'  \item{DI}{SpatRaster, stars object or data frame. Dissimilarity index of newdata}
#'  \item{LPD}{SpatRaster, stars object or data frame. Local Point Density of newdata.}
#'  \item{AOA}{SpatRaster, stars object or data frame. Area of Applicability of newdata. AOA has values 0 (outside AOA) and 1 (inside AOA)}
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
#' library(terra)
#' library(caret)
#' library(viridis)
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
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea, model, LPD = TRUE, maxLPD = 1)
#' plot(AOA)
#'
#' ####
#' #The AOA can also be calculated without a trained model.
#' #All variables are weighted equally in this case:
#' ####
#' AOA <- aoa(studyArea,train=trainDat,variables=variables)
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
#' prediction <- predict(studyArea,lrn$model,na.rm=TRUE)
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
                useWeight=TRUE,
                LPD = FALSE,
                maxLPD = 1) {

  # handling of different raster formats
  as_stars <- FALSE
  leading_digit <- any(grepl("^{1}[0-9]",names(newdata)))

  if (inherits(newdata, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    newdata <- methods::as(newdata, "SpatRaster")
    as_stars <- TRUE
  }
  if (inherits(newdata, "Raster")) {
   # if (!requireNamespace("raster", quietly = TRUE))
  #    stop("package raster required: install that first")
    message("Raster will soon not longer be supported. Use terra or stars instead")
    newdata <- methods::as(newdata, "SpatRaster")
  }

  calc_LPD <- LPD
  # validate maxLPD input
  if (LPD == TRUE) {
    if (is.numeric(maxLPD)) {
      if (maxLPD <= 0) {
        stop("maxLPD can not be negative or equal to 0. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
      }
      if (maxLPD <= 1) {
        if (inherits(model, "train")) {
          maxLPD <- round(maxLPD * as.integer(length(model$trainingData[[1]])))
        } else if (!is.null(train)) {
          maxLPD <- round(maxLPD * as.integer(length(train[[1]])))
        }
        if (maxLPD <= 1) {
          stop("The percentage you provided for maxLPD is too small.")
        }
      }
      if (maxLPD > 1) {
        if (maxLPD %% 1 == 0) {
          maxLPD <- as.integer(maxLPD)
        } else if (maxLPD %% 1 != 0) {
          stop("If maxLPD is bigger than 0, it should be a whole number. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
        }
      }
      if ((maxLPD > length(if (inherits(model, "train")) { model$trainingData[[1]] } else if (!is.null(train)) { train[[1]] })) || maxLPD %% 1 != 0) {
        stop("maxLPD can not be bigger than the number of training samples. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
      }
    } else {
      stop("maxLPD must be a number. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
    }
  }


  # if not provided, compute train DI
  if(!inherits(trainDI, "trainDI")) {
    message("No trainDI provided.")
    trainDI <- trainDI(model, train, variables, weight, CVtest, CVtrain, method, useWeight, LPD)
  }

  if (calc_LPD == TRUE) {
    # maxLPD <- trainDI$avrgLPD
    trainDI$maxLPD <- maxLPD
  }


  # check if variables are in newdata
  if(any(trainDI$variables %in% names(newdata)==FALSE)){
    if(leading_digit){
      stop("names of newdata start with leading digits, automatically added 'X' results in mismatching names of train data in the model")
    }
    stop("names of newdata don't match names of train data in the model")
  }


  # Prepare output as either as RasterLayer or vector:
  out <- NA
  if (inherits(newdata, "SpatRaster")){
    out <- newdata[[1]]
    names(out) <- "DI"
  }



  #### order data:
  if (inherits(newdata, "SpatRaster")){
    if (any(is.factor(newdata))){
      newdata[[which(is.factor(newdata))]] <- as.numeric(newdata[[which(is.factor(newdata))]])
    }
    newdata <- terra::as.data.frame(newdata,na.rm=FALSE)
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
  okrows <- which(apply(newdata, 1, function(x)
    all(!is.na(x))))
  newdataCC <- newdata[okrows, ]

  if (method == "MD") {
    if (dim(train_scaled)[2] == 1) {
      S <- matrix(stats::var(train_scaled), 1, 1)
      newdataCC <- as.matrix(newdataCC, ncol = 1)
    } else {
      S <- stats::cov(train_scaled)
    }
    S_inv <- MASS::ginv(S)
  }

  if (calc_LPD == FALSE) {
    message("Computing DI of new data...")
    mindist <- rep(NA, nrow(newdata))
    mindist[okrows] <-
      .mindistfun(newdataCC, train_scaled, method, S_inv)
    DI_out <- mindist / trainDI$trainDist_avrgmean
  }

  # if (calc_LPD == TRUE) {
  #   message("Computing LPD of newdata...")
  #
  #   knndist <- matrix(NA, nrow(newdata), maxLPD)
  #   knndist[okrows,] <- .knndistfun(newdataCC, train_scaled, method, S_inv, maxLPD = maxLPD)
  #
  #   DI_out_knndist <- knndist / trainDI$trainDist_avrgmean
  #   DI_out <- c(DI_out_knndist[,1])
  #
  #   # start_time <- Sys.time()
  #   LPD_out <-
  #     apply(DI_out_knndist, 1, function(row)
  #       sum(row < trainDI$threshold))
  #   # end_time <- Sys.time()
  #
  #   # set maxLPD to max of LPD_out if
  #   realMaxLPD <- max(LPD_out, na.rm = T)
  #   if (maxLPD > realMaxLPD) {
  #     if (inherits(maxLPD, c("numeric", "integer"))) {
  #       message("Your specified maxLPD is bigger than the real maxLPD of you predictor data.")
  #     }
  #     message(paste("maxLPD is set to", realMaxLPD))
  #     trainDI$maxLPD <- realMaxLPD
  #   }
  # }

  if (calc_LPD == TRUE) {
    message("Computing DI and LPD of new data...")

    pb <- txtProgressBar(min = 0,
                         max = nrow(newdataCC),
                         style = 3)

    DI_out <- rep(NA, nrow(newdata))
    LPD_out <- rep(NA, nrow(newdata))
    for (i in seq(nrow(newdataCC))) {
      knnDist  <- .knndistfun(t(matrix(newdataCC[i,])), train_scaled, method, S_inv, maxLPD = maxLPD)
      knnDI <- knnDist / trainDI$trainDist_avrgmean
      knnDI <- c(knnDI)

      DI_out[okrows[i]] <- knnDI[1]
      LPD_out[okrows[i]] <- sum(knnDI < trainDI$threshold)
      setTxtProgressBar(pb, i)
    }

    close(pb)

    # set maxLPD to max of LPD_out if
    realMaxLPD <- max(LPD_out, na.rm = T)
    if (maxLPD > realMaxLPD) {
      if (inherits(maxLPD, c("numeric", "integer"))) {
        message("Your specified maxLPD is bigger than the real maxLPD of you predictor data.")
      }
      message(paste("maxLPD is set to", realMaxLPD))
      trainDI$maxLPD <- realMaxLPD
    }
  }


  message("Computing AOA...")

  #### Create Mask for AOA and return statistics
  if (inherits(out, "SpatRaster")) {
    terra::values(out) <- DI_out

    AOA <- out
    terra::values(AOA) <- 1
    AOA[out > trainDI$thres] <- 0
    AOA <- terra::mask(AOA, out)
    names(AOA) = "AOA"

    if (calc_LPD == TRUE) {
      LPD <- out
      terra::values(LPD) <- LPD_out
      names(LPD) = "LPD"
    }


    # handling of different raster formats.
    if (as_stars) {
      out <- stars::st_as_stars(out)
      AOA <- stars::st_as_stars(AOA)

      if (calc_LPD == TRUE) {
        LPD <- stars::st_as_stars(LPD)
      }
    }

  } else{
    out <- DI_out
    AOA <- rep(1, length(out))
    AOA[out > trainDI$thres] <- 0

    if (calc_LPD == TRUE) {
      LPD <- LPD_out
    }
  }


  # used in old versions of the AOA. eventually remove the attributes
#  attributes(AOA)$aoa_stats <- list("Mean_train" = trainDI$trainDist_avrgmean,
#                                    "threshold" = trainDI$thres)
#  attributes(AOA)$TrainDI <- trainDI$trainDI

  if (calc_LPD == FALSE) {
    result <- list(
      parameters = trainDI,
      DI = out,
      AOA = AOA
    )
  } else {
    result <- list(
      parameters = trainDI,
      DI = out,
      AOA = AOA,
      LPD = LPD
    )
  }

  message("Finished!")

  class(result) <- "aoa"
  return(result)
}


.knndistfun <-
  function (point,
            reference,
            method,
            S_inv = NULL,
            maxLPD = maxLPD) {
    if (method == "L2") {
      # Euclidean Distance
      return(FNN::knnx.dist(reference, point, k = maxLPD))
    } else if (method == "MD") {
      return(t(sapply(1:dim(point)[1],
                      function(y)
                        sort(sapply(1:dim(reference)[1],
                                    function(x)
                                      sqrt(t(point[y, ] - reference[x, ]) %*% S_inv %*% (point[y, ] - reference[x,]) )))[1:maxLPD])))
    }
  }

