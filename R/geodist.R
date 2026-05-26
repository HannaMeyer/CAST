#' Calculate euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Calculates nearest neighbor distances in geographic space or feature space between training data as well as between prediction locations and training data.
#' Optional, the nearest neighbor distances between test data and training data or between different CV folds is computed.
#' @param x object of class sf, training data locations
#' @param modeldomain SpatRaster, stars or sf object defining the prediction area (see Details)
#' @param dist_space "geographical", "feature" or "time". Should the distance be computed in geographic space, in the normalized multivariate predictor space or in temporal space? (see Details)
#' @param CVtest optional. list or vector. #' @param cvfolds optional. list or vector. Either a list with the length of the number of cross-validation folds
#' where each element contains the row indices of the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm
#' @param CVtrain optional. A list, where each element contains the data points used for training during the cross validation iteration.
#' Only required if CVtrain is not the opposite of CVtest. Relevant if some data points are excluded, e.g. when using \code{\link{nndm}}.
#' @param testdata optional. object of class sf: Point data used for independent validation. May already include the predictor values if `dist_space`=feature.
#' @param preddata optional. object of class sf: Point data indicating the locations within the modeldomain to be used as target prediction points. Useful when the prediction objective is a subset of
#' locations within the modeldomain rather than the whole area. May already include the predictor values if `dist_space`="feature".
#' @param samplesize numeric. How many prediction samples should be used?
#' @param sampling character. How to draw prediction samples? See \link[sf]{st_sample} for modeldomains that are sf objects and \link[terra]{spatSample} for raster objects.
#' Use sampling = "Fibonacci" for global applications (raster objects will be transformed to polygons in this case).
#' @param variables character vector defining the predictor variables used if dist_space="feature". If not provided all variables included in modeldomain are used.
#' @param time_var optional. character. Column that indicates the date. Only used if dist_space="time".
#' @param time_unit optional. Character. Unit for temporal distances See ?difftime.Only used if dist_space="time".
#' @param dist_fun character. Automatically detected if `dist_space`="geographical". For geographical (long/lat) coordinates,
#' `dist_fun` is set to "great_circle", while "euclidean" distances are used for projected coordinats.
#' For `dist_space`="feature", `dist_fun` currently covers `euclidean` (default), `gower` and `mahalanobis`.
#' `mahalanobis` takes into account correlation between predictor values. While `euclidean` and `mahalanobis` only work with numerical variables,
#' `gower` also works with mixed data including numerical and categorical variables.
#' For `dist_space`="time", currently only the absolute difference (`abs_time`) is implemented.
#' @param scale_vars boolean. Should variables be scaled? Only for `dist_space`="feature". 
#' Calculating Gower distances already includes scaling, and manually rescale the data is redundant. 
#' For other distances (Mahalanobis, Euclidean), scaling the data is important. Thus, TRUE by default.
#' @param cvtrain deprecated. Use `CVtrain` instead.
#' @param cvfolds deprecated. Use `CVtest` instead.
#' @param type deprecated. Use `dist_space` instead.
#' @param timevar deprecated. Use `time_var` instead.
#' @return A data.frame containing the distances. Unit of returned geographic distances is meters. attributes contain W statistic between prediction area and either sample data, CV folds or test data. See details.
#' @details The modeldomain is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by samplesize) from the spatial extent (if no `preddata` are supplied).
#'     If `dist_space` = "feature", the argument modeldomain has to be a raster and include predictors. The only exception is when the provided training data and preddata already include the predictor values.
#'     If not provided they are extracted from the modeldomain raster. If some predictors are categorical (i.e., of class factor or character), gower distances will be used.
#'     W statistic describes the match between the distributions. See Linnenbrink et al (2024) for further details.
#' @note See Meyer and Pebesma (2022) for an application of this plotting function
#' @seealso \code{\link{nndm}} \code{\link{knndm}}
#' @import ggplot2
#' @author Hanna Meyer, Edzer Pebesma, Marvin Ludwig, Jan Linnenbrink
#' @examples
#' \dontrun{
#' library(CAST)
#' library(sf)
#' library(terra)
#' library(caret)
#' library(rnaturalearth)
#' library(ggplot2)
#'
#' data(splotdata)
#' studyArea <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf")
#'
#' ########### Distance between training data and new data:
#' dist <- geodist(splotdata, studyArea)
#' # With density functions
#' plot(dist)
#' # Or ECDFs (relevant for nndm and knnmd methods)
#' plot(dist, stat="ecdf")
#'
#' ########### Distance between training data, new data and test data (here Chile):
#' plot(splotdata[,"Country"])
#' dist <- geodist(splotdata[splotdata$Country != "Chile",], studyArea,
#'                 testdata = splotdata[splotdata$Country == "Chile",])
#' plot(dist)
#'
#' ########### Distance between training data, new data and CV folds:
#' folds <- createFolds(1:nrow(splotdata), k=3, returnTrain=FALSE)
#' dist <- geodist(x=splotdata, modeldomain=studyArea, CVtest=folds)
#' # Using density functions
#' plot(dist)
#' # Using ECDFs (relevant for nndm and knnmd methods)
#' plot(dist, stat="ecdf")
#'
#' ########### Distances in the feature space:
#' predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))
#' dist <- geodist(x = splotdata,
#'                 modeldomain = predictors,
#'                 dist_space = "feature",
#'                 variables = c("bio_1","bio_12", "elev"))
#' plot(dist)
#'
#' dist <- geodist(x = splotdata[splotdata$Country != "Chile",],
#'                 modeldomain = predictors,
#'                 testdata = splotdata[splotdata$Country == "Chile",],
#'                 dist_space = "feature",
#'                 variables=c("bio_1","bio_12", "elev"))
#' plot(dist)
#'
#'############Distances in temporal space
#' library(lubridate)
#' library(ggplot2)
#' data(cookfarm)
#' dat <- st_as_sf(cookfarm,coords=c("Easting","Northing"))
#' st_crs(dat) <- 26911
#' trainDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2010,]
#' predictionDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2011,]
#' trainDat$week <- lubridate::week(trainDat$Date)
#' CVtest <- CreateSpacetimeFolds(trainDat,time_var = "week")
#'
#' dist <- geodist(trainDat,preddata = predictionDat,CVtest = CVtest$indexOut,
#'    dist_space="time",time_unit="days")
#' plot(dist)+ xlim(0,10)
#'
#'
#' ############ Example for a random global dataset
#' ############ (refer to figure in Meyer and Pebesma 2022)
#'
#' ### Define prediction area (here: global):
#' ee <- st_crs("+proj=eqearth")
#' co <- ne_countries(returnclass = "sf")
#' co.ee <- st_transform(co, ee)
#'
#' ### Simulate a spatial random sample
#' ### (alternatively replace pts_random by a real sampling dataset (see Meyer and Pebesma 2022):
#' sf_use_s2(FALSE)
#' pts_random <- st_sample(co.ee, 2000, exact=FALSE)
#'
#' ### See points on the map:
#' ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
#'   geom_sf(data = pts_random, color = "#F8766D",size=0.5, shape=3) +
#'   guides(fill = "none", col = "none") +
#'   labs(x = NULL, y = NULL)
#'
#' ### plot distances:
#' dist <- geodist(pts_random,co.ee)
#' plot(dist) + scale_x_log10(labels=round)
#'
#'
#'
#'
#'}
#' @export
geodist <- function(
  x,
  modeldomain = NULL,
  dist_space = "geographical",
  CVtest = NULL,
  CVtrain = NULL,
  testdata = NULL,
  preddata = NULL,
  samplesize = 2000,
  sampling = "regular",
  variables = NULL,
  time_var = NULL,
  time_unit = "auto",
  dist_fun = "euclidean",
  scale_vars = TRUE,
  cvtrain = NULL,
  cvfolds = NULL,
  type = NULL,
  timevar = NULL
){

  ## Guard against old parameter names / deprecated parameters
  if (!is.null(cvtrain)) {
    warning("Argument 'cvtrain' is deprecated. Please use 'CVtrain' instead.",
            call. = FALSE)
    CVtrain <- cvtrain
  }

  if (!is.null(cvfolds)) {
    warning("Argument 'cvfolds' is deprecated. Please use 'CVtest' instead.",
            call. = FALSE)
    CVtest <- cvfolds
  }

  if (!is.null(type)) {
    warning("Argument 'type' is deprecated. Please use 'dist_space' instead.",
            call. = FALSE)
    dist_space <- type
  }

  if (!is.null(timevar)) {
    warning("Argument 'timevar' is deprecated. Please use 'time_var' instead.",
            call. = FALSE)
    time_var <- timevar
  }

  if (dist_space == "geo") dist_space <- "geographical"

  ## Check that dist_space was correctly defined
  if (!dist_space %in% c("geographical", "feature", "time")) {
    stop("dist_space must be one of 'geographical', 'feature' or 'time'")
  }
  if (!(dist_fun %in% c("euclidean", "mahalanobis", "gower", "great_circle", "abs_time"))) {
    stop("dist_fun must be one of 'euclidean', 'mahalanobis', 'gower' or 'great_circle'")
  }
  if(dist_space == "time" && dist_fun != "abs_time"){
    warning("Temporal space only supports 'abs_time' distances.")
    dist_fun <- "abs_time"
  } 
  if(dist_space == "feature" && dist_fun == "great_circle") {
    stop("Great-circle distances only work with in geographical space.")
  }
  if(dist_space == "geographical" && dist_fun %in% c("mahalanobis", "gower")) {
    stop("Mahalanobis and Gower distances only work in feature space.")
  }

  ## CVtrain
  if(!is.null(CVtrain) && !is.list(CVtrain)) {
    stop("CVtrain has to be a list of indices")
  }
  if(!is.null(CVtest) && !is.list(CVtest)) {
    stop("CVtest has to be a list of indices")
  }
  if(!is.null(CVtrain) && is.null(CVtest)) {
    message("CVtest was inferred as the opposite of CVtrain")
    n_flds <- max(unlist(CVtrain))
    CVtest <- lapply(CVtrain, function(train_idx) setdiff(seq_len(n_flds), train_idx))
  }
  if (!is.null(CVtest) && is.null(CVtrain)) {
    n_flds <- max(unlist(CVtest))
    CVtrain <- lapply(CVtest, function(test_idx) setdiff(seq_len(n_flds), test_idx))
  }
  folds <- NULL
  if (!is.null(CVtest) || !is.null(CVtrain)) folds <- list(train = CVtrain, test = CVtest)

  ## Check for different raster formats
  if (inherits(modeldomain, "Raster")) {
    modeldomain <- methods::as(modeldomain,"SpatRaster")
  }
  if (inherits(modeldomain, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE)) {
      stop("package stars required: install that first")
    }
    modeldomain <- methods::as(modeldomain, "SpatRaster")
  }

  # Retrieve variable names if not provided
  variables <- .get_ref_vars(variables, x, modeldomain, preddata)

  # Sample prediction points from the study area if not supplied
  if (is.null(preddata)) {
    pred_points <- sampleFromArea(
      modeldomain = modeldomain, 
      samplesize = samplesize, 
      dist_space = dist_space, 
      variables = variables, 
      sampling = sampling)
  } else {
    pred_points <- preddata
  }

  # Retrieve the reference CRS of preddata (or modeldomain)
  ref_crs <- .get_ref_crs(pred_points, modeldomain, x)
  # transform training points to the reference crs if necessary
  if (sf::st_crs(x) != ref_crs) {
    message("Transforming training data CRS to match the preddata/modeldomain")
    x <- sf::st_transform(x, ref_crs)
  }
  # transform test points to the reference crs if necessary
  if (!is.null(testdata) && sf::st_crs(testdata) != sf::st_crs(x)) {
    message("Transforming test data CRS to match the preddata/modeldomain")
    testdata <- sf::st_transform(testdata, sf::st_crs(x))
  }

  # now hand off to the appropriate distance processor function
  if (dist_space == "geographical") {
    message("Calculating distances in geographic space...")
    nnds <- .geo_processor(
      x = x, 
      pred_points = pred_points,
      testdata = testdata, 
      dist_fun = dist_fun,
      folds = folds
    )
  } else if (dist_space == "time") {
    message("Calculating distances in temporal space...")
    nnds <- .time_processor(
      x = x,
      pred_points = pred_points,
      testdata = testdata,
      dist_fun = dist_fun,
      folds = folds,
      time_var = time_var,
      time_unit = time_unit
    )
  } else if (dist_space == "feature") {
    message("Calculating distances in feature space...")
    nnds <- .feature_processor(
      x = x, 
      pred_points = pred_points, 
      testdata = testdata,
      dist_fun = dist_fun,
      folds = folds, 
      variables = variables,
      scale_vars = scale_vars,
      modeldomain = modeldomain
    )
  }
  # Post-processing ----------

  class(nnds) <- c("geodist", class(nnds))
  attr(nnds, "dist_space") <- dist_space
  if (dist_space == "time") attr(nnds, "unit") <- time_unit

  attr(nnds, "W_sample") <- twosamples::wass_stat(
    nnds[nnds$what == "sample-to-sample", "dist"],
    nnds[nnds$what == "prediction-to-sample", "dist"]
  )

  if (!is.null(testdata)) {
    attr(nnds, "W_test") <- twosamples::wass_stat(
      nnds[nnds$what == "test-to-sample", "dist"],
      nnds[nnds$what == "prediction-to-sample", "dist"]
    )
  }

  if (!is.null(CVtest)) {
    attr(nnds, "W_CV") <- twosamples::wass_stat(
      nnds[nnds$what == "CV-distances", "dist"],
      nnds[nnds$what == "prediction-to-sample", "dist"]
    )
  }
  nnds
}

.get_ref_crs <- function(preddata, modeldomain, x) {
  crs <- sf::st_crs(x)
  if (!is.null(modeldomain)) {
    crs <- sf::st_crs(modeldomain)
  }
  if (!is.null(preddata)) {
    crs <- sf::st_crs(preddata)
  }
  return(crs)
}


# Processor functions for different space types ----------
.geo_processor <- function(
  x, 
  pred_points = NULL, 
  testdata = NULL,
  dist_fun = c("euclidean", "great_circle"),
  folds = NULL) {
  
  dist_fun <- match.arg(dist_fun)
  
  islonglat <- if (is.na(sf::st_crs(x))) {
    warning("Missing CRS of the modeldomain or prediction points. Assuming projected CRS.")
    FALSE
  } else {
    sf::st_is_longlat(sf::st_crs(x))
  }

  if (islonglat) {
    message("Calculating great-circle distances in geographic space (longlat coordinates).")
    dist_fun <- "great_circle"
  } else {
    message("Calculating euclidean distances in geographic space (projected coordinates).")
    dist_fun <- "euclidean"
  }

  # Calculate NNDs between training points
  dists <- .compute_nnds(
    x = x, 
    pred_points = pred_points, 
    testdata = testdata, 
    folds = folds,
    .dist = .geo_dist, 
    args = list(dist_fun = dist_fun)
  )
  dists$dist_type <- "geographical"
  return(dists)
}

.time_processor <- function(
  x, 
  pred_points = NULL, 
  testdata = NULL,
  dist_fun = "abs_time", 
  folds = NULL,
  time_var = NULL, 
  time_unit = "auto") {
  
  dist_fun <- match.arg(dist_fun)

  if (is.null(time_var)) {
    time_var <- names(which(vapply(x, lubridate::is.Date, logical(1))))
    if (length(time_var) == 0) {
      stop("No time variable found. Please specify the time variable using the 'time_var' argument.")
    } else if (length(time_var) > 1) {
      warning("Multiple time variables found. Using the first one: ", time_var[1])
      time_var <- time_var[1]
    }
    message("time variable selected: ", time_var)
  }

  dists <- .compute_nnds(
    x = x, 
    pred_points = pred_points, 
    testdata = testdata, 
    folds = folds, 
    .dist = .time_dist, 
    args = list(time_var = time_var, time_unit = time_unit)
  )
  dists$dist_type <- "time"
  return(dists)
}

.feature_processor <- function(
  x, 
  pred_points = NULL,
  testdata = NULL,
  dist_fun = c("euclidean", "mahalanobis", "gower"),
  folds = NULL,
  variables = NULL,
  scale_vars = TRUE,
  modeldomain = NULL) {
  
  
  dist_fun <- match.arg(dist_fun)
  if (is.null(pred_points) && !inherits(modeldomain, "SpatRaster")) {
    stop("For feature space, modeldomain must be a SpatRaster or preddata must be supplied.")
  }
  # retrive variable names
  variables <- .get_ref_vars(variables, x, modeldomain, pred_points)
  # Extract predictor values from the modeldomain if they are not attached
  x <- .extract_predictors(x, modeldomain, variables)
  pred_points <- .extract_predictors(pred_points, modeldomain, variables)
  if (!is.null(testdata)) {
    testdata <- .extract_predictors(testdata, modeldomain, variables)
  }

  # Drop the geometry of the training and prediction (and test) points
  x <- sf::st_drop_geometry(x)
  pred_points <- sf::st_drop_geometry(pred_points)
  if (!is.null(testdata)) testdata <- sf::st_drop_geometry(testdata)

  # Detect categorical variables
  cat_cols <- .get_categorical_variables(x, variables)
  num_cols <- setdiff(variables, cat_cols)

  if (length(cat_cols) > 0 && dist_fun != "gower") {
      stop("Only 'gower' distances are allowed for categorical variables.")
  }

  if (isTRUE(scale_vars)) {
    if (length(num_cols) > 0) { # we only need to scale if there are numerical variables
      scale_attr <- attributes(scale(x[, num_cols, drop=FALSE]))
      center <- scale_attr$`scaled:center`
      scale <- scale_attr$`scaled:scale`
      x[, num_cols] <- scale(x[, num_cols, drop=FALSE])
      pred_points[, num_cols] <- scale(pred_points[, num_cols, drop=FALSE], center = center, scale = scale)
      if (!is.null(testdata)) {
        testdata[, num_cols] <- scale(testdata[, num_cols, drop=FALSE], center = center, scale = scale)
      }
    }
  }

  dists <- .compute_nnds(
    x = x, 
    pred_points = pred_points, 
    testdata = testdata, 
    folds = folds, 
    .dist = .feat_dist, 
    args = list(dist_fun = dist_fun)
  )
  dists$dist_type <- "feature"
  return(dists)
}

.get_ref_vars <- function(vars, x, modeldomain, preddata) {
  if (!is.null(vars)) return(vars)
  if (!is.null(modeldomain)) return(names(modeldomain))
  if (!is.null(preddata)) return(names(preddata))
  return(names(x))
}

.extract_predictors <- function(x, modeldomain, variables) {
  if (any(!variables %in% names(x))) {
    message("Extracting predictors from modeldomain")
    x <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(x), bind = TRUE, na.rm = TRUE))
  }
  x <- x[, variables]
  return(x)
}

##### atomical distance functions ----------
.geo_dist <- function(
  x, 
  y = NULL, 
  dist_fun = NULL) {
  
  if(dist_fun == "great_circle") {
    dist_mat <- if(is.null(y)) sf::st_distance(x) else sf::st_distance(y, x)
    units(dist_mat) <- NULL
    if(is.null(y)) diag(dist_mat) <- NA
    min_d <- apply(dist_mat, 1, min, na.rm=TRUE)
  } 

  if (dist_fun == "euclidean") {
    coords_x <- sf::st_coordinates(x)[ ,1:2]
    coords_y <- if(!is.null(y)) sf::st_coordinates(y)[ ,1:2] else NULL
    min_d <- .knndist(query = coords_y, reference = coords_x, k = 1, dist_fun = dist_fun)
  }

  return(min_d)
}

.time_dist <- function(
  x, 
  y = NULL,
  time_var, 
  time_unit = "auto"){
  
  time_x <- sf::st_drop_geometry(x)[ ,time_var]
  time_y <- if(!is.null(y)) sf::st_drop_geometry(y)[ ,time_var] else time_x
  
  if (time_unit == "auto") {
    time_unit <- units(difftime(time_x, time_x))
  }

  dist_mat <- abs(outer(time_y, time_x, FUN = function(a,b) as.numeric(difftime(a,b,units=time_unit))))
  if(is.null(y)) diag(dist_mat) <- NA
  min_d <- apply(dist_mat, 1, min, na.rm = TRUE)

  return(min_d)
}

.feat_dist <- function(
  x, 
  y = NULL,
  dist_fun) {
  min_d <- .knndist(
    query = y, 
    reference = x, 
    k = 1, 
    dist_fun = dist_fun, 
    offset = 0 # TODO: check if this is really correct: do we need to ignore self distances when y is NULL?
  )
  return(min_d)
}


# Helper function: retrieve near neighbors distances ----------
.compute_nnds <- function(
  x, 
  pred_points, 
  testdata, 
  folds = NULL,
  .dist,
  args
) {
  # sample to sample distances
  s2s <- data.frame(dist = do.call(.dist, c(list(x = x, y = NULL), args)))
  s2s$what <- "sample-to-sample"
  # prediction to sample distances
  p2s <- data.frame(dist = do.call(.dist, c(list(x = x, y = pred_points), args)))
  p2s$what <- "prediction-to-sample"
  # target to sample distances (if testdata is provided)
  t2s <- NULL
  if (!is.null(testdata)) {
    t2s <- data.frame(dist = do.call(.dist, c(list(x = x, y = testdata), args)))
    t2s$what <- "test-to-sample"
  }
  # CV fold distances (if CVtest is provided)
  cvdist <- NULL
  if (!is.null(folds) && !is.null(folds$test)) {
    cvdist <- data.frame(dist = .cv_dist(x = x, folds = folds, .dist = .dist, args = args))
    cvdist$what <- "CV-distances"
  }
  # Combine distances into a single data frame
  dists <- rbind(s2s, p2s)
  if (!is.null(t2s)) dists <- rbind(dists, t2s)
  if (!is.null(cvdist)) dists <- rbind(dists, cvdist)
  return(dists)
}

# Helper function: call .dist either with or without CV folds ---------
.cv_dist <- function(x, folds, .dist, args) {
  # Get total number of observations
  n <- nrow(x)
  # Initialize result vector to preserve position ordering
  all_cv_dists <- rep(NA_real_, n)
  # Process each fold and assign distances to correct positions
  for (i in seq_along(folds$test)) {
    test_idx <- folds$test[[i]]
    train_idx <- folds$train[[i]]
    fold_dists <- do.call(.dist, c(list(x = x[train_idx, , drop=FALSE], y = x[test_idx, , drop=FALSE]), args))
    # Assign to correct positions in result vector
    all_cv_dists[test_idx] <- fold_dists
  }
  return(all_cv_dists)
}

## Helper function: sample prediction points from the prediction area ----------
sampleFromArea <- function(modeldomain, samplesize, dist_space, variables, sampling){

  # Sample points from a Raster
  if(inherits(modeldomain, "Raster")){
    modeldomain <- terra::rast(modeldomain)
  }

  if(inherits(modeldomain, "SpatRaster") && sampling != "Fibonacci") {

    if(samplesize>terra::ncell(modeldomain)){
      samplesize <- terra::ncell(modeldomain)
      message(paste0("samplesize for new data shouldn't be larger than number of pixels.
              Samplesize was reduced to ",terra::ncell(modeldomain)))
    }
    #create mask to sample from:
    template <- modeldomain[[1]]
    template <- terra::classify(template, cbind(-Inf, Inf, 1), right=FALSE)
    # draw samples using terra
    message(paste0("Sampling ", samplesize, " prediction locations from the modeldomain raster."))
    predictionloc <- terra::spatSample(template, size = samplesize, method = sampling, as.points = TRUE, na.rm = TRUE, values = FALSE) |> 
      sf::st_as_sf()

  } else if(inherits(modeldomain, "SpatRaster") && sampling == "Fibonacci") {
    message("Converting raster to polygon for Fibonacci sampling")
    # Use first layer as template
    raster_template <- modeldomain[[1]]
    # Convert non-NA cells to polygons
    template <- sf::st_as_sf(terra::as.polygons(raster_template, na.rm = TRUE))
    # Sample points from the polygon
    message(paste0("Sampling ", samplesize, " prediction locations from the modeldomain vector."))
    predictionloc <- suppressMessages(sf::st_sample(template, size = samplesize, type = sampling)) |> 
      sf::st_set_crs(sf::st_crs(modeldomain))

  } else {
    # Sample points from a Polygon
    message(paste0("Sampling ", samplesize, " prediction locations from the modeldomain vector."))
    predictionloc <- suppressMessages(sf::st_sample(modeldomain, size = samplesize, type = sampling)) |> 
      sf::st_set_crs(sf::st_crs(modeldomain))
  }
  return(predictionloc)

}
