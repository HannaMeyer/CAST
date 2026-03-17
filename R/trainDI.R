#' Calculate Dissimilarity Index of training data
#' @description
#' This function estimates the Dissimilarity Index (DI)
#' within the training data set used for a prediction model.
#' Optionally, the local point density can also be calculated.
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
#' @param CVtest list or vector. Either a list with the length of the number of cross-validation folds where each element contains the row indices of the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point.
#' Only required if no model is given.
#' @param CVtrain list. Each element contains the data points used for training during the cross validation iteration (i.e. held back data).
#' Only required if no model is given and only required if CVtrain is not the opposite of CVtest (i.e. if a data point is not used for testing, it is used for training).
#' Relevant if some data points are excluded, e.g. when using \code{\link{nndm}}.
#' @param dist_fun Character. Method used for distance calculation. Currently, euclidean and mahalanobis distance are implemented but only euclidean is tested.
#' @param useWeight Logical. Only if a model is given. Weight variables according to importance in the model?
#' @param useCV Logical. Only if a model is given. Use the CV folds to calculate the DI threshold?
#' @param LPD Logical. Indicates whether the local point density should be calculated or not.
#' @param chunk_size Integer. Number of training points to be processed in each chunk when calculating distances. Decreasing this number can help to reduce memory usage but increases runtime.
#' @param verbose Logical. Print progress or not?
#' @param method Deprecated. Use dist_fun instead.
#' @param algorithm Deprecated. Use dist_fun instead.
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
#' Hanna Meyer, Marvin Ludwig, Fabian Schumacher
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
#' library(CAST)
#'
#' # prepare sample data:
#' data("splotdata")
#' splotdata = st_drop_geometry(splotdata)
#'
#' # train a model:
#' set.seed(100)
#' model <- caret::train(splotdata[,6:16],
#'                       splotdata$Species_richness,
#'                       importance=TRUE, tuneLength=1, ntree = 15, method = "rf",
#'                       trControl = trainControl(method="cv", number=5, savePredictions=T))
#' # variable importance is used for scaling predictors
#' plot(varImp(model,scale=FALSE))
#'
#' # calculate the DI of the trained model:
#' DI = trainDI(model=model)
#' plot(DI)
#'
#' #...or calculate the DI and LPD of the trained model:
#' # DI = trainDI(model=model, LPD = TRUE)
#'
#' # the DI can now be used to compute the AOA (here with LPD):
#' studyArea = rast(system.file("extdata/predictors_chile.tif", package = "CAST"))
#' AOA = aoa(studyArea, model = model, trainDI = DI, LPD = TRUE, maxLPD = 1)
#' print(AOA)
#' plot(AOA)
#' plot(AOA$AOA)
#' plot(AOA$LPD)
#' }
#'


trainDI <- function(model = NULL,
                    train = NULL,
                    variables = "all",
                    weight = NULL,
                    CVtest = NULL,
                    CVtrain = NULL,
                    dist_fun=c("euclidean", "mahalanobis", "gower"),
                    useWeight = TRUE,
                    useCV =TRUE,
                    LPD = FALSE,
                    chunk_size = 1000L,
                    verbose = TRUE,
                    method,
                    algorithm) {
  
  if (!missing(method) || !missing(algorithm)) {
    warning("The 'method' and 'algorithm' parameters are deprecated. Please use 'dist_fun' instead.")
  }
  # backwards compatibility if user code specifies model or weight as NA instead of NULL
  if (!is.null(model) && is.na(model)[1]) model <- NULL
  if (!is.null(weight) && is.na(weight)[1]) weight <- NULL

  dist_fun <- match.arg(dist_fun)
  # get parameters if they are not provided in function call-----
  if(is.null(train)){train = aoa_get_train(model)}
  if(length(variables) == 1 && variables == "all"){
      variables = aoa_get_variables(variables, model, train)
  }
  if(is.null(weight)){
    if(useWeight){
      weight = aoa_get_weights(model, variables = variables)
    }else{
      message("variable are not weighted. see ?aoa")
      weight <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
      names(weight) <- variables
    }
  }else{
    weight <- user_weights(weight, variables)
  }

  # get CV folds from model or from parameters
  folds <-  aoa_get_folds(model,CVtrain,CVtest,useCV)
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
  scaleparam <- attributes(train)

  # make sure all variables have variance
  if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
    stop("some variables in train seem to have no variance")
  }

  # multiply train data with variable weights (from variable importance)
  if(!inherits(weight, "error")&!is.null(unlist(weight))){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }

  # calculate average mean distance between training data
  train_dists <- chunked_dist(train=train, CVtrain=CVtrain, CVtest=CVtest, 
    dist_fun=dist_fun, chunk_size = chunk_size, verbose = verbose)

  # Dissimilarity Index defined as the minimum distance to the nearest neighbour 
  # divided by the average distance to all other points
  trainDist_avrgmean =  mean(train_dists$trainDist_avrg, na.rm = TRUE)
  TrainDI <- train_dists$trainDist_min / trainDist_avrgmean
  thres <- aoa_threshold(TrainDI)

  # prepare aoa return object
  aoa_results = list(
    train = train_backup,
    weight = weight,
    variables = variables,
    catvars = catupdate$catvars,
    scaleparam = scaleparam,
    trainDist_avrg = train_dists$trainDist_avrg,
    trainDist_avrgmean = trainDist_avrgmean,
    trainDI = TrainDI,
    threshold = thres,
    method = dist_fun
  )

  # calculate trainLPD and avrgLPD according to the CV folds of specified
  if (LPD) {
    trainLPD <- chunked_lpd(train=train, CVtrain=CVtrain, CVtest=CVtest, 
        dist_fun=dist_fun, train_mean = trainDist_avrgmean, 
        threshold = thres, chunk_size = chunk_size, verbose = verbose) 
    aoa_results$trainLPD <- trainLPD
    aoa_results$avrgLPD <- round(mean(trainLPD))
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
    weight <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
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
  if(sum(weight)==0){
    stop("all weights are <=0, hence no variable is used. Check variable importance of the model, define weights manually or set useWeight=FALSE")
  }
  return(weight)

}



# check user weight input
# make sure this function outputs a data.frame with
# one row and columns named after the variables

user_weights = function(weight, variables){

  # list input support
  if(inherits(weight, "list")){
    # check if all list entries are in variables
    weight = as.data.frame(weight)
  }


  #check if manually given weights are correct. otherwise ignore (set to 1):
  if(nrow(weight)!=1  || !all(variables %in% names(weight))){
    message("variable weights are not correctly specified and will be ignored. See ?aoa")
    weight <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
    names(weight) <- variables
  }
  weight <- weight[,na.omit(match(variables, names(weight)))]
  if (any(weight<0)){
    weight[weight<0]<-0
    message("negative weights were set to 0")
  }

  return(weight)

}




# Get trainingdata from train object

aoa_get_train <- function(model){

  train <- as.data.frame(model$trainingData)
  return(train)


}


# Get folds from train object


aoa_get_folds <- function(model, CVtrain, CVtest, useCV){
  ### if folds are to be extracted from the model:
  if (useCV&!is.null(model)){
    if(tolower(model$control$method)!="cv"){
      message("note: Either no model was given or no CV was used for model training. The DI threshold is therefore based on all training data")
    }else{
      CVtest <- model$control$indexOut
      CVtrain <- model$control$index
    }
  }
  ### if folds are specified manually:
  if(is.null(model)){

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
  if(!is.null(model)&useCV==FALSE){
    message("note: useCV is set to FALSE. The DI threshold is therefore based on all training data")
    CVtrain <- NULL
    CVtest <- NULL
  }
  return(list(CVtrain,CVtest))
}


# Get variables from train object

aoa_get_variables <- function(variables, model, train){
  if (length(variables) > 1) {
    return(variables)
  }
  if(variables == "all"){
    if(is.null(model)[1]){
      variables <- names(train)
    }else{
      variables <- names(model$trainingData)[-which(names(model$trainingData)==".outcome")]
    }
  }
  return(variables)
}

# helper to derive DI threshold based on the distribution of DI in the training data
aoa_threshold <- function(trainDI){
  threshold_quantile <- stats::quantile(trainDI, 0.75, na.rm=TRUE)
  threshold_iqr <- (1.5 * stats::IQR(trainDI, na.rm=T))
  thres <- threshold_quantile + threshold_iqr
  # make sure that the threshold is not larger than the maximum DI in the training data
  thres <- min(thres, max(trainDI, na.rm=T)) 
  return(thres)
}

# small helper for distances and lpd calculation in chunks to reduce memory usage
chunked_dist <- function(
  train, 
  CVtrain, 
  CVtest, 
  dist_fun, 
  chunk_size = 1000L, 
  verbose = TRUE) {
  res <- chunked_apply(
    train = train,
    CVtrain = CVtrain,
    CVtest = CVtest,
    dist_fun = dist_fun,
    chunk_size = chunk_size,
    verbose = verbose,
    calc_fun = calc_dist
  )
  # list with min distance, average distance, index of nearest neighbour for each observation
  return(res)
}

chunked_lpd <- function(
  train, 
  CVtrain, 
  CVtest, 
  dist_fun, 
  train_mean, 
  threshold, 
  chunk_size = 1000L, 
  verbose = TRUE) {
  lpd <- chunked_apply(
    train = train,
    CVtrain = CVtrain,
    CVtest = CVtest,
    dist_fun = dist_fun,
    chunk_size = chunk_size,
    verbose = verbose,
    calc_fun = calc_lpd,
    train_mean = train_mean,
    threshold = threshold
  )
  # vector of local point density (counts of neighbours below DI threshold) per observation
  return(lpd)
}


#' Calculate chunked KNN distances for training data
#'
#' Internal helper that computes per-query-row average distance to all reference
#' points and the minimum distance to the nearest allowed training neighbour.
#' Distances are computed with knndist and masked via mask_dist_mat to remove
#' self- and within-fold neighbours.
#'
#' @param reference numeric matrix or data.frame of reference (scaled & weighted) rows.
#' @param query numeric matrix or data.frame of query rows (subset of reference).
#' @param CVtrain list of training indices per CV fold or NULL.
#' @param CVtest list of testing indices per CV fold or NULL.
#' @param dist_fun character; distance method forwarded to knndist.
#' @param ids integer vector mapping rows of query to row indices in reference.
#'
#' @return A named list with
#'   \item{trainDist_min}{numeric vector of minimum distances (per query row).}
#'   \item{trainDist_avrg}{numeric vector of average distances (per query row).}
#'   \item{trainDist_indices}{integer vector of index of nearest neighbour (per query row).}
#'
#' @keywords internal
#' @noRd
calc_dist <- function(
  reference, 
  query = NULL, 
  CVtrain, 
  CVtest, 
  dist_fun, 
  ids){
  
  dist_mat <- knndist(query = query, reference = reference,
                      k = nrow(reference) - 1, dist_fun = dist_fun,
                      return_distmat = TRUE)

  # first, we only mask the observations themselves and retrieve the average distance to all other points
  dist_mat <- mask_dist_mat(dist_mat = dist_mat, ids = ids)
  avg <- apply(dist_mat, 1, mean, na.rm = TRUE)
  # now we mask within-fold observations and retrieve the minimum distance to the closest point
  dist_mat <- mask_dist_mat(dist_mat = dist_mat, ids = ids, CVtest = CVtest, CVtrain = CVtrain)
  # TODO: consider retrieving th K-th nearest neighbour and its distance
  dist <- apply(dist_mat, 1, min, na.rm = TRUE)
  idx <- apply(dist_mat, 1, function(x) { which(x == min(x, na.rm = TRUE))[1] })

  list(
    # vector of minimum distances to nearest neighbour (per training point)
    trainDist_min = dist,
    # vector of average distances to all other points (per training point)
    trainDist_avrg = avg,
    # vector of indices of nearest neighbour (per training point)
    trainDist_indices = idx)
}


#' Calculate local point density (LPD) for a chunk
#'
#' Internal helper that computes, for each query row, the number of reference
#' points whose (masked) DI = distance / train_mean is below a given threshold.
#' Distances are computed with knndist and masked via mask_dist_mat.
#'
#' @param reference numeric matrix or data.frame of reference (scaled & weighted) rows.
#' @param query numeric matrix or data.frame of query rows (subset of reference).
#' @param CVtrain list of training indices per CV fold or NULL.
#' @param CVtest list of testing indices per CV fold or NULL.
#' @param dist_fun character; distance method forwarded to knndist.
#' @param train_mean numeric scalar; normalization factor (mean training distance).
#' @param threshold numeric scalar; DI threshold used to count neighbours.
#' @param ids integer vector mapping rows of query to row indices in reference.
#'
#' @return Integer vector with the local point density (counts) for each query row.
#'
#' @keywords internal
#' @noRd
calc_lpd <- function(
  reference, 
  query = NULL, 
  CVtrain, 
  CVtest, 
  dist_fun, 
  train_mean, 
  threshold, 
  ids){
  
  dist_mat <- knndist(query = query, reference = reference,
                      k = nrow(reference) - 1, dist_fun = dist_fun,
                      return_distmat = TRUE)
  # mask within-fold observations and self-distances to get the local point density
  dist_mat <- mask_dist_mat(dist_mat = dist_mat, ids = ids, 
    CVtest = CVtest, CVtrain = CVtrain)
  # convert distances to DI by normalizing with the global mean distance, 
  # then count how many neighbours are below the DI threshold
  di_mat <- dist_mat / train_mean
  trainLPD <- as.integer(rowSums(di_mat < threshold, na.rm = TRUE))
  trainLPD
}

mask_dist_mat <- function(dist_mat, ids, CVtest = NULL, CVtrain = NULL) {
  if (is.null(dim(dist_mat))) return(dist_mat) # nothing to do for vectors
  # dist_mat here is an NxM matrix of distances from N query points (rows) to M reference points (columns)

  # Set self-distances to NA: each row corresponds to ids[j] in the reference
  for (j in seq_along(ids)) {
    row_idx <- j
    col_idx <- ids[j]
    dist_mat[row_idx, col_idx] <- NA
  }

  # Mask within-fold distances when CV folds are provided
  # first -> retrieve the testing fold for each observation
  # second -> keep only the training samples of the respective observation 
  # and set everything else to NA
  if (!is.null(CVtrain) && !is.null(CVtest)) {
    for (j in seq_along(ids)) {
      sample_idx <- ids[j]
      # find the testing fold that contains the sample
      whichfold <- which(vapply(CVtest, function(x) any(x == sample_idx), logical(1))) 
      if(length(whichfold) > 1L) { # raise if a sample is used for testing more than once
        stop("a datapoint is used for testing in more than one fold. currently this option is not implemented")
      }
      if (length(whichfold) == 0L) { # if never used in testing, we ignore the sample completely
        dist_mat[j, ] <- NA
      } else { # otherwise, we set all non-training samples to NA
        train_ids <- CVtrain[[whichfold]]
        if (length(train_ids) > 0) {
          not_train_ids <- setdiff(seq_len(ncol(dist_mat)), train_ids)
          dist_mat[j, not_train_ids] <- NA
        }
      }
    }
  }
  dist_mat
}


#' Apply a calculation function over training data in row-wise chunks
#'
#' Internal utility that splits `train` into chunks of size `chunk_size`,
#' calls `calc_fun(reference = train, query = chunk, ..., ids = chunk_ids)`,
#' and combines results. If `calc_fun` returns a named list for each chunk,
#' results are concatenated per name and returned as a named list; otherwise
#' chunk outputs are concatenated into a single vector.
#'
#' @param train data.frame or matrix of training data (rows = samples).
#' @param CVtrain list of training indices per CV fold or NULL.
#' @param CVtest list of testing indices per CV fold or NULL.
#' @param dist_fun character; distance method passed through to calc_fun.
#' @param chunk_size integer; number of rows per chunk.
#' @param verbose logical; print progress messages.
#' @param calc_fun function to call for each chunk. Must accept at least
#'   arguments (reference, query, CVtrain, CVtest, dist_fun, ids, ...).
#' @param ... additional arguments forwarded to calc_fun.
#'
#' @return Either a named list (if calc_fun returns named lists per chunk)
#'   with each element concatenated across chunks, or a concatenated vector
#'   of chunk results.
#'
#' @keywords internal
#' @noRd
chunked_apply <- function(
  train,
  CVtrain,
  CVtest,
  dist_fun,
  chunk_size = 1000L,
  verbose = TRUE,
  calc_fun,
  ...
) {
  # split train into chunks of size chunk_size and apply calc_fun to each chunk
  n <- nrow(train)
  if (n == 0) return(list())
  n_chunks <- ceiling(n / chunk_size)
  chunk_ids <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  results <- vector("list", n_chunks)

  # TODO: consider parallelization of the chunk processing,
  # e.g. with future.apply::future_lapply or furry::future_map
  for (i in seq_len(n_chunks)) {
    if (verbose) message(sprintf("Processing chunk %d of %d", i, n_chunks))
    # ensure row-subset keeps columns; pass ids for masking inside calc_fun
    results[[i]] <- calc_fun(
      reference = train,
      query     = train[chunk_ids[[i]], , drop = FALSE],
      CVtrain   = CVtrain,
      CVtest    = CVtest,
      dist_fun  = dist_fun,
      ids       = chunk_ids[[i]],
      ...
    )
  }

  first <- results[[1]]

  # If each chunk returns a named list -> combine by name
  if (is.list(first) && !is.null(names(first))) {
    out <-vector("list", length(first))
    names(out) <- names(first)
    for (nm in names(first)) {
      out[[nm]] <- unlist(lapply(results, function(x) x[[nm]]))
    }
    return(out)
  }

  # Otherwise assume each chunk returned a vector -> concatenate
  unlist(results)
}
