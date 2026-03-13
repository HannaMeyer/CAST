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
#' Optionally, the local point density is calculated which indicates the number of similar training data points up to the DI threshold.
#' @param newdata A SpatRaster, stars object or data.frame containing the data
#' the model was meant to make predictions for.
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds.
#' See examples for the case that no model is available or for models trained via e.g. mlr3.
#' @param trainDI A trainDI object. Optional if \code{\link{trainDI}} was calculated beforehand.
#' @param train A data.frame containing the data used for model training. Optional. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Optional. Only required if no model is given.
#' @param ... ignored
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the model are used or if no model is given then of the train dataset.
#' @param CVtest list or vector. Either a list with the length of the number of cross-validation folds where each element contains the row indices of the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point.
#' Only required if no model is given.
#' @param CVtrain list. Each element contains the data points used for training during the cross validation iteration (i.e. held back data).
#' Only required if no model is given and only required if CVtrain is not the opposite of CVtest (i.e. if a data point is not used for testing, it is used for training).
#' Relevant if some data points are excluded, e.g. when using \code{\link{nndm}}.
#' @param method Character. Method used for distance calculation. Currently euclidean distance (L2) and Mahalanobis distance (MD) are implemented but only L2 is tested. Note that MD takes considerably longer.
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param useCV Logical. Only if a model is given. Use the CV folds to calculate the DI threshold?
#' @param LPD Logical. Indicates whether the local point density should be calculated or not.
#' @param maxLPD numeric or integer. Only if \code{LPD = TRUE}. Number of nearest neighbors to be considered for the calculation of the LPD. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples. CAUTION! If not all training samples are considered, a fitted relationship between LPD and error metric will not make sense (@seealso \code{\link{DItoErrormetric}})
#' @param indices logical. Calculate indices of the training data points that are responsible for the LPD of a new prediction location? Output is a matrix with the dimensions num(raster_cells) x maxLPD. Each row holds the indices of the training data points that are relevant for the specific LPD value at that location. Can be used in combination with exploreAOA(aoa) function from the \href{https://github.com/fab-scm/CASTvis}{CASTvis package} for a better visual interpretation of the results. Note that the matrix can be quite big for examples with a high resolution and a larger number of training samples, which can cause memory issues.
#' @param parallel Logical. Parallelization the process. Only possible if LPD = TRUE. Can reduce computation time significantly.
#' @param cores Integer or Character. Number of cores to use for the the parallelization. You can use "auto" to set your cores to \code{detectCores()/2} (see \code{\link[parallel]{detectCores}}).
#' @param verbose Logical. Print progress or not?
#' @details The Dissimilarity Index (DI), the Local Data Point Density (LPD) and the corresponding Area of Applicability (AOA) are calculated.
#' If variables are factors, dummy variables are created prior to weighting and distance calculation.
#'
#' Interpretation of results: If a location is very similar to the properties
#' of the training data it will have a low distance in the predictor variable space
#' (DI towards 0) while locations that are very different in their properties
#' will have a high DI. For easier interpretation see \code{\link{normalize_DI}}
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
#' @importFrom parallel detectCores makeForkCluster clusterExport parLapply stopCluster
#'
#' @author
#' Hanna Meyer, Fabian Schumacher
#' @references Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' Methods in Ecology and Evolution 12: 1620-1633. \doi{10.1111/2041-210X.13650}
#'
#' Schumacher, F., Knoth, C., Ludwig, M., Meyer, H. (2025):
#' Estimation of local training data point densities to support the assessment
#' of spatial prediction uncertainty. Geosci. Model Dev., 18, 10185–10202. \doi{10.5194/gmd-18-10185-2025}
#'
#' @seealso \code{\link{trainDI}}, \code{\link{normalize_DI}}, \code{\link{errorProfiles}}
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(caret)
#'
#' # prepare sample data:
#' data(cookfarm)
#' dat <- aggregate(cookfarm[,c("VW","Easting","Northing")],
#'    by=list(as.character(cookfarm$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"),crs=26911)
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
#' AOA <- aoa(studyArea, model)
#' plot(AOA)
#' plot(AOA$AOA)
#' #... or if preferred calculate the aoa and the LPD of the study area:
#' AOA <- aoa(studyArea, model, LPD = TRUE)
#' plot(AOA$LPD)
#'
#' #note that it is not required to use Random Forests. The method is model agnostic.
#' # Let's chnage to SVM:
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW, method="svmRadial", importance=TRUE, tuneLength=1,
#' trControl=trainControl(method="cv",number=5,savePredictions=T))
#' AOA <- aoa(studyArea, model, LPD = TRUE)
#' plot(AOA$LPD)
#'
#' ####
#' #The AOA can also be calculated without a trained model.
#' #All variables are weighted equally in this case:
#' ####
#'
#' AOA <- aoa(studyArea,train=trainDat,variables=variables)
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

#' @export
#' @name aoa
aoa = function(newdata, model=NA, ...) UseMethod("aoa") # the generic

#' @export
#' @name aoa
aoa.stars = function(newdata, model=NA, ...) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    res = aoa(methods::as(newdata, "SpatRaster"), model=model, ...)
	# convert back to stars...
    res$DI <- stars::st_as_stars(res$DI)
    res$AOA <- stars::st_as_stars(res$AOA)
    if (!is.null(res$LPD)) {
        res$LPD <- stars::st_as_stars(res$LPD)
    }
	res
}

#' @export
#' @name aoa
aoa.Raster = function(newdata, model=NA, ...) {
    stop("Raster will soon not longer be supported. Use terra or stars instead")
    aoa(methods::as(newdata, "SpatRaster"), model=model, ...)
}

#' @export
#' @name aoa
aoa.SpatRaster = function(newdata, model=NA, ...) {
  #### order data:
  out_template = newdata[[1]]
  if (any(is.factor(newdata)))
    newdata[[which(is.factor(newdata))]] <- as.numeric(newdata[[which(is.factor(newdata))]])
  # call the data.frame method:
  res = aoa(terra::as.data.frame(newdata, na.rm = FALSE), model=model, ...)

  # convert DI:
  DI = out_template
  terra::values(DI) <- res$DI
  DI <- terra::mask(DI, out_template)
  names(DI) = "DI"
  res$DI = DI

  # convert AOA:
  AOA = out_template
  terra::values(AOA) <- res$AOA
  AOA <- terra::mask(AOA, out_template)
  names(AOA) = "AOA"
  res$AOA = AOA

  if (!is.null(res$LPD)) {
    # convert LPD:
    LPD <- out_template
    terra::values(LPD) <- res$LPD
    names(LPD) = "LPD"
	res$LPD = LPD
  }
  res
}


#' @export
#' @name aoa
aoa.data.frame <- function(newdata,
                model=NA,
				...,
                trainDI = NA,
                train=NULL,
                weight=NA,
                variables="all",
                CVtest=NULL,
                CVtrain=NULL,
                method="euclidean",
                useWeight=TRUE,
                useCV=TRUE,
                LPD = FALSE,
                maxLPD = 1,
                indices = FALSE,
                parallel = FALSE,
                cores = 4,
                verbose = TRUE) {

  leading_digit <- any(grepl("^{1}[0-9]",names(newdata)))

  calc_LPD <- LPD
  # validate maxLPD input
  if (LPD == TRUE) {
    if (is.numeric(maxLPD)) {
      if (maxLPD <= 0) {
        stop("maxLPD cannot be negative or equal to 0. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
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
        stop("maxLPD cannot be bigger than the number of training samples. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
      }
    } else {
      stop("maxLPD must be a number. Either define a number between 0 and 1 to use a percentage of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
    }
  }

  if (parallel & Sys.info()["sysname"] != "Linux") {
    stop("Paralellization only works for UNIX-alike systems. Please use single core computation.")
  }

  # if not provided, compute train DI
  if(!inherits(trainDI, "trainDI")) {
    if (verbose) {
      message("No trainDI provided.")
    }
    trainDI <- trainDI(
      model=model, 
      train=train, 
      variables=variables, 
      weight=weight,
      CVtest=CVtest, 
      CVtrain=CVtrain, 
      method=method, 
      useWeight=useWeight, 
      useCV=useCV, 
      LPD=LPD, 
      verbose=verbose)
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


  # TODO:
  # Prepare output as either as RasterLayer or vector:
  out <- NA
  if (inherits(newdata, "SpatRaster")){
    out <- newdata[[1]]
    names(out) <- "DI"
  }

  newdata <- newdata[,na.omit(match(trainDI$variables, names(newdata))),drop = FALSE]

  ## Handling of categorical predictors:
  catvars <- trainDI$catvars
  if (!inherits(catvars,"error")&length(catvars)>0){
    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
	  # x[,catvar] is a vector if x is a data.frame, but a tibble if x is a tibble
	  # x[[catvar]] extracts the catvar variable for both data.frame and tibble
      trainDI$train[[catvar]]<-droplevels(trainDI$train[[catvar]])
      newdata[[catvar]] <- factor(newdata[[catvar]])
      newdata[!newdata[[catvar]]%in%unique(trainDI$train[[catvar]]),catvar] <- NA
      newdata[[catvar]] <- droplevels(newdata[[catvar]])
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
    tmpnames <- names(newdata)#!!!!!
    newdata <- sapply(1:ncol(newdata),function(x){
      newdata[,x]*unlist(trainDI$weight[x])
    })
    names(newdata)<-tmpnames#!!!!
  }


  # rescale and reweight train data
  train_scaled <- scale(trainDI$train,
                        center = trainDI$scaleparam$`scaled:center`,
                        scale = trainDI$scaleparam$`scaled:scale`)

  train_scaled <- sapply(1:ncol(train_scaled),function(x){train_scaled[,x]*unlist(trainDI$weight[x])})


  # Distance Calculation ---------
  okrows <- which(rowSums(is.na(newdata)) == 0)
  newdataCC <- newdata[okrows, ,drop=F]

  if (!calc_LPD) {
    if (verbose) {
      message("Computing DI of new data...")
    }
    mindist <- rep(NA, nrow(newdata))
    mindist[okrows] <-
      .knndistfun(query=newdataCC, reference=train_scaled, k=1, dist_fun=method)
    DI_out <- mindist / trainDI$trainDist_avrgmean
  }

  if (calc_LPD) {
    if (verbose) {
      message("Computing DI and LPD of new data...")
    }

    DI_out <- rep(NA, nrow(newdata))
    LPD_out <- rep(NA, nrow(newdata))
    if (indices) {
        Indices_out <- matrix(NA, nrow = nrow(newdataCC), ncol = maxLPD)
    }

    knnDist  <- .knndistfun(query=newdataCC, reference=train_scaled, k=maxLPD, dist_fun=method)
    knnDI <- knnDist / trainDI$trainDist_avrgmean
    DI_out[okrows] <- knnDI[ ,1]
    LPD_out[okrows] <- rowSums(knnDI < trainDI$threshold)
    if (indices) {
      if (LPD_out[okrows] > 0) {
        knnIndex  <- .knndistfun(query=newdataCC, reference=train_scaled, k = LPD_out[okrows[i]], dist_fun=method, distance = FALSE)
        Indices_out[i,1:LPD_out[okrows]] <- as.numeric(knnIndex)
      }
    }

    # set maxLPD to max of LPD_out if
    realMaxLPD <- max(LPD_out, na.rm = T)
    if (maxLPD > realMaxLPD) {
      if (inherits(maxLPD, c("numeric", "integer")) && verbose) {
        message("Your specified maxLPD is bigger than the real maxLPD of you predictor data.")
      }
      if (verbose) {
        message(paste("maxLPD is set to", realMaxLPD))
      }
      trainDI$maxLPD <- realMaxLPD
    }

    if (indices) {
      Indices_out <- Indices_out[,1:trainDI$maxLPD]
      rownames(Indices_out) <- okrows
    }
  }

  if (verbose) {
    message("Computing AOA...")
  }

  result <- list(
    parameters = trainDI,
    DI = DI_out,
    AOA = ifelse(DI_out > trainDI$threshold, 0L, 1L)
  )

  if (calc_LPD == TRUE) {
    result$LPD <- LPD_out
    if (indices) {
      result$indices <- Indices_out
    }
  }

  if (verbose) {
    message("Finished!")
  }

  class(result) <- "aoa"
  return(result)
}

.process_row <- function(row) {
  knnDist <- .knndistfun(query=row, reference=train_scaled, k=maxLPD, dist_fun=method)
  knnDI <- knnDist / trainDIdat$trainDist_avrgmean
  knnDI <- c(knnDI)

  DI_out_i <- knnDI[1]
  LPD_out_i <- sum(knnDI < trainDIdat$threshold)

  # return if indices not to be calculated
  return(list(DI_out_i = DI_out_i,
              LPD_out_i = LPD_out_i
  ))
}

# Tell R CMD check these variables are fine
utils::globalVariables(
  c(
    "train_scaled",
    "method",
    "trainDIdat",
    "maxLPD",
    "indices"
    )
  )


