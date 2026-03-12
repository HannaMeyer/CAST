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
#' @param algorithm see \code{\link[FNN]{knnx.dist}} and \code{\link[FNN]{knnx.index}}
#' @param dist_fun character. Currently covers `euclidean` (default), `gower`, `mahalanobis`, `great_circle` and `abs_time`.
#' `gower` and `mahalanobis` only work with `dist_space`="feature", while `great_circle` only works with `dist_space`="geographical". 
#' `mahalanobis` takes into account correlation between predictor values. While `euclidean` and `mahalanobis` only work with numerical variables, 
#' `gower` also works with mixed data including numerical and categorical variables.
#' For `dist_space`="time", currently only the absolute difference (`abs_time`) is implemented.
#' For the geographical space, `great_circle` covers lon/lat coordinates, whereas `euclidean` only works with projected coordinates.
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
  algorithm = "brute",
  dist_fun = "euclidean",
  scale_vars = TRUE,
  cvtrain = NULL,
  cvfolds = NULL,
  type = NULL,
  timevar = NULL
){

  # 1. Input validation & normalization ----------

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
  if(dist_space == "feature" && dist_fun == "great_circle") stop("Great-circle distances only work with in geographical space.")
  if(dist_space == "geographical" && dist_fun %in% c("mahalanobis", "gower")) stop("Mahalanobis and Gower distances only work in feature space.")

  ## CVtrain
  if(!is.null(CVtrain) && !is.list(CVtrain)) stop("CVtrain has to be a list of indices")
  if(is.null(CVtest) && !is.null(CVtrain)) {
    message("CVtest was inferred as the opposite of CVtrain")
    n_flds <- max(unlist(CVtrain))
    # Generate CVtest as the complement of CVtrain for each fold
    CVtest <- lapply(CVtrain, function(train_idx) setdiff(seq_len(n_flds), train_idx))
  }


  ## Check for different raster formats
  if (inherits(modeldomain, "Raster")) {
    modeldomain <- methods::as(modeldomain,"SpatRaster")
  }
  if (inherits(modeldomain, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    modeldomain <- methods::as(modeldomain, "SpatRaster")
  }

  ## Check for time variable and unit
  if (dist_space == "time") {
    if (is.null(time_var)) {
      time_var <- names(which(vapply(x, lubridate::is.Date, logical(1))))
      message("time variable selected: ", time_var)
    }
    if (time_unit == "auto") {
      time_unit <- units(difftime(
        sf::st_drop_geometry(x)[, time_var],
        sf::st_drop_geometry(x)[, time_var]
      ))
    }
  }

  # 2. CRS harmonization ----------

  # Retrieve the CRS of preddata (or modeldomain) and use it as a reference
  ref_crs <- if (!is.null(preddata)) {
    sf::st_crs(preddata)
  } else if (!is.null(modeldomain) && inherits(modeldomain, "sf")) {
    sf::st_crs(modeldomain)
  } else if (!is.null(modeldomain) && inherits(modeldomain, "SpatRaster")) {
    sf::st_crs(modeldomain)
  } else {
    sf::st_crs(x)
  }

  # transform training points to the CRS of preddata (or modeldomain)
  if (!is.na(ref_crs) && !is.na(sf::st_crs(x)) && sf::st_crs(x) != ref_crs) {
    message("Transforming training data CRS to match the preddata/modeldomain")
    x <- sf::st_transform(x, ref_crs)
  }

  # transform test points to the CRS of preddata (or modeldomain)
  if (!is.null(testdata) &&
      !is.na(sf::st_crs(testdata)) &&
      !is.na(sf::st_crs(x)) &&
      sf::st_crs(testdata) != sf::st_crs(x)) {
    message("Transforming test data CRS to match the preddata/modeldomain")
    testdata <- sf::st_transform(testdata, sf::st_crs(x))
  }

  # Check if coordinates are longitutude/latitude
  islonglat <- if (is.na(ref_crs)) {
    warning("Missing CRS of the modeldomain or prediction points. Assuming projected CRS.")
    FALSE
  } else {
    sf::st_is_longlat(ref_crs)
  }

  # 3. Prediction point generation ----------

  # Sample prediction points from the study area (only if no preddata are supplied)
  pred_points <- if (is.null(preddata)) {
    sampleFromArea(modeldomain, samplesize, dist_space, variables, sampling)
  } else {
    preddata
  }


  # 4. Feature-space preparation (if needed) ----------

  catVars <- NULL

  if (dist_space == "feature") {

    if (is.null(preddata) && !inherits(modeldomain, "SpatRaster")) {
      stop("For feature space, modeldomain must be a SpatRaster or preddata must be supplied.")
    }

    # If no variable names are given, retrieve them from the preddata or modeldomain (the latter is preferred)
    if (is.null(variables) && !is.null(modeldomain)) {
       variables <- names(modeldomain)
    } else if (is.null(variables) && !is.null(preddata) && is.null(modeldomain)) {
       variables <- names(preddata)
    }

    # Extract predictor values from the modeldomain if they are not attached
    if (any(!variables %in% names(x))) {
      message("Extracting predictors from modeldomain")
      x <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(x), bind = TRUE, na.rm = TRUE))
    }
    x <- x[, variables]

    if (!is.null(testdata) && any(!variables %in% names(testdata))) {
      testdata <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(testdata), bind = TRUE, na.rm = TRUE))
    }
    if (!is.null(testdata)) testdata <- testdata[, variables]

    if (any(!variables %in% names(pred_points))) {
      pred_points <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(pred_points), bind = TRUE, na.rm = TRUE))
    }
    pred_points <- pred_points[, variables]

    # Detect categorical variables
    catVars <- names(x)[vapply(x, function(z) inherits(z, c("factor", "character")), logical(1))]
    if (length(catVars) == 0) catVars <- NULL
    # Drop the geometry of the training and prediction (and test) points
    x <- sf::st_drop_geometry(x)
    pred_points <- sf::st_drop_geometry(pred_points)
    if (!is.null(testdata)) testdata <- sf::st_drop_geometry(testdata)

    # Optionally scale the training, prediction and test points
    if(isTRUE(scale_vars)) {
            scaleparam <- attributes(scale(x))
            x <- data.frame(scale(x))
            pred_points <- data.frame(scale(pred_points,
                                            center = scaleparam$`scaled:center`,
                                            scale = scaleparam$`scaled:scale`))
            if (!is.null(testdata)) {
              testdata <- data.frame(scale(testdata,
                                          center = scaleparam$`scaled:center`,
                                          scale = scaleparam$`scaled:scale`))
            }
    }
    
  }

  # 5. Distance computation ----------

  # Some more input checks
  if(dist_space == "feature") {
    if (!is.null(catVars) && dist_fun != "gower") {
        stop("Only 'gower' distances are allowed for categorical variables.")
    }
  } else if(dist_space == "geographical") {
    if(isTRUE(islonglat) && dist_fun != "great_circle")  {
      stop("Only 'great_circle' distances are allowed for lon/lat coordinates.")
    }
  }


  # Calculate NNDs between training points
  s2s <- compute_NND(
    x=x, y=NULL, dist_space = dist_space, dist_type_label = "sample-to-sample",  dist_fun = dist_fun,
    CVtest = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
  )

  # Calculate NNDs between prediction and training points
  p2s <- compute_NND(
    x=x, y=pred_points, dist_space = dist_space, dist_type_label = "prediction-to-sample",  dist_fun = dist_fun,
    CVtest = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
  )
  dists <- rbind(s2s, p2s)

  # Calculate NNDs between test points and training points
  if (!is.null(testdata)) {
    t2s <- compute_NND(
      x=x, y=testdata, dist_space = dist_space, dist_type_label = "test-to-sample",  dist_fun = dist_fun,
      CVtest = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
    )
    dists <- rbind(dists, t2s)
  }
  # Calculate NNDs between CV folds
  if (!is.null(CVtest)) {
    cvdist <- compute_NND(
      x=x, y=NULL, dist_space = dist_space, dist_type_label = "CV-distances",  dist_fun = dist_fun,
      CVtest = CVtest, CVtrain = CVtrain, time_var = time_var, time_unit = time_unit, algorithm = algorithm
    )
    dists <- rbind(dists, cvdist)
  }

  # 6. Post-processing ----------

  class(dists) <- c("geodist", class(dists))
  attr(dists, "dist_space") <- dist_space
  if (dist_space == "time") attr(dists, "unit") <- time_unit

  attr(dists, "W_sample") <- twosamples::wass_stat(
    dists[dists$what == "sample-to-sample", "dist"],
    dists[dists$what == "prediction-to-sample", "dist"]
  )

  if (!is.null(testdata)) {
    attr(dists, "W_test") <- twosamples::wass_stat(
      dists[dists$what == "test-to-sample", "dist"],
      dists[dists$what == "prediction-to-sample", "dist"]
    )
  }

  if (!is.null(CVtest)) {
    attr(dists, "W_CV") <- twosamples::wass_stat(
      dists[dists$what == "CV-distances", "dist"],
      dists[dists$what == "prediction-to-sample", "dist"]
    )
  }

  dists
}


## Distance calculation based on space ----------
compute_NND <- function(x, y = NULL, dist_space = c("geographical","feature","time"), 
                      dist_type_label = "sample-to-sample", 
                      dist_fun = c("euclidean","mahalanobis","gower","great_circle", "abs_time"), 
                      CVtest = NULL, CVtrain = NULL, time_var = NULL, time_unit = "auto", algorithm = "brute") {
  
  dist_space <- match.arg(dist_space)
  dist_fun <- match.arg(dist_fun)

  # Stop when encountering invalid parameter combinations
  if(!is.null(y) && !is.null(CVtest)) stop("Currently, `CVtest` and `y` cannot be specified simultaniously.")

  # Calculate nearest neighbor distances in geographical space
  if(dist_space == "feature"){
    # Calculate nearest neighbor distances in feature space
    reference <- sf::st_drop_geometry(x)
    query <- if(!is.null(y)) sf::st_drop_geometry(y) else NULL

    if(is.null(CVtest)) {
      min_d <- .knndistfun(query = query, reference = reference, k = 1, dist_fun = dist_fun, distance = TRUE)
      } else {
      min_d <- cv_distances(reference, CVtest = CVtest, CVtrain = CVtrain, dist_fun = dist_fun)
      }
 
  } else if(dist_space == "geographical") {
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
    # Extract coordinates from x (and y)
    coords_x <- sf::st_coordinates(x)[,1:2]
    coords_y <- if(!is.null(y)) sf::st_coordinates(y)[,1:2] else NULL

    if(dist_fun == "great_circle") {
      if(is.null(CVtest)) {
        distmat <- if(is.null(y)) sf::st_distance(x) else sf::st_distance(y, x)
        units(distmat) <- NULL
        if(is.null(y)) diag(distmat) <- NA
        min_d <- apply(distmat, 1, min, na.rm=TRUE)
      } else {
        min_d <- cv_distances(x, CVtest = CVtest, CVtrain = CVtrain, dist_fun = dist_fun)
      }
    } else { # euclidean
        if(is.null(CVtest)) {
          min_d <- .knndistfun(query = coords_y, reference = coords_x, k = 1, dist_fun = dist_fun, distance = TRUE)
        } else {
          min_d <- cv_distances(coords_x, CVtest = CVtest, CVtrain = CVtrain, dist_fun = dist_fun)
        }
    }
  } else if(dist_space == "time"){
    # Calculate nearest neighbor distances in temporal space
    if(dist_fun == "abs_time") {
      time_x <- sf::st_drop_geometry(x)[,time_var]
      if(is.null(CVtest)) {
        time_y <- if(!is.null(y)) sf::st_drop_geometry(y)[,time_var] else time_x
        dmat <- abs(outer(time_y, time_x, FUN = function(a,b) as.numeric(difftime(a,b,units=time_unit))))
        if(is.null(y)) diag(dmat) <- Inf
        min_d <- apply(dmat, 1, min)
      } else {
        min_d <- cv_distances(time_x, CVtest = CVtest, time_unit = time_unit, CVtrain = CVtrain, dist_fun = dist_fun)
      }
    }
  }
  
  data.frame(dist = as.vector(min_d),
             what = factor(dist_type_label),
             dist_type = dist_space)
}


## Helper function: Compute out-of-fold NN distance ----------
cv_distances <- function(x, CVtest, CVtrain = NULL, dist_fun = "euclidean", time_unit = "auto") {

  # define length of NND vector
  n <- if(inherits(x, c("Date","POSIXct","POSIXt"))) length(x) else nrow(x)
  alldist <- rep(NA, n)

  # Convert CVtest and CVtrain to vectors (knndm supplies CVtest in vector format, so we can't only support lists of indices)
  if(inherits(CVtest, "list")) {
    CVtest_v <- rep(NA_integer_, n)
    for(k in seq_along(CVtest)) CVtest_v[CVtest[[k]]] <- k
  } else CVtest_v <- CVtest
    
  # Calculate distances between CV folds
  for(f in unique(CVtest_v)){
    
    test_idx <- which(CVtest_v == f)
    
    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(CVtest_v != f)
    }

    if(inherits(x, c("Date","POSIXct","POSIXt"))) {
        tr_train <- x[train_idx]
        tr_test  <- x[test_idx]
      } else {
        tr_train <- x[train_idx,, drop=FALSE]
        tr_test  <- x[test_idx,, drop=FALSE]
      }

    if(dist_fun %in% c("euclidean", "mahalanobis", "gower")) {
      alldist[test_idx] <- .knndistfun(query = tr_test, reference = tr_train, k = 1, dist_fun = dist_fun, distance = TRUE)
    } else if(dist_fun == "great_circle") {
        distmat <- sf::st_distance(x)
        units(distmat) <- NULL
        diag(distmat) <- NA 
        alldist[test_idx] <- apply(distmat[test_idx, train_idx, drop=FALSE], 1, min)
    } else if(dist_fun == "abs_time") {
      diffs <- outer(tr_test, tr_train, function(x,y) abs(as.numeric(difftime(x, y, units=time_unit))))
      alldist[test_idx] <- apply(diffs, 1, min)
    }    
  }
  alldist
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

.knndistfun <- function(
  reference, 
  query = NULL,
  k = 1, 
  dist_fun = c("euclidean", "mahalanobis", "gower"), 
  distance = TRUE) {
  
  dist_fun <- match.arg(dist_fun)

  if (dist_fun == "gower") {
    # requires scaling of the numerical variables, while categorical variables are not scaled.
    num_vars <- vapply(reference, is.numeric, logical(1))
    ref_num <- reference[, num_vars, drop = FALSE]
    mins <- sapply(ref_num, min, na.rm = TRUE)
    maxs <- sapply(ref_num, max, na.rm = TRUE)
    range <- maxs - mins

    reference[, num_vars] <- as.data.frame(
        sweep(sweep(reference[, num_vars, drop = FALSE], 2, mins, "-"), 2, range, "/"),
        stringsAsFactors = FALSE
    )
    
    if (!is.null(query)) {
      query[, num_vars] <- as.data.frame(
        sweep(sweep(query[, num_vars, drop = FALSE], 2, mins, "-"), 2, range, "/"),
        stringsAsFactors = FALSE
      )
    }
    query <- .numeric_fct(query)
    reference <- .numeric_fct(reference)
  }
  
  if (inherits(query, "numeric")) {
    query <- matrix(query, nrow = 1)
  }
  if (inherits(query, "data.frame")) {
    query <- as.matrix(query)
  }
  if (inherits(reference, "numeric")) {
    reference <- matrix(reference, nrow = 1)
  }
  if (inherits(reference, "data.frame")) {
    reference <- as.matrix(reference)
  }

  if (dist_fun == "mahalanobis") {
    # For Mahalanobis distance, we need to compute the inverse covariance matrix 
    # we then transform that reference and query to calculate the L2 distance in 
    # the transformed space, which is equivalent to the Mahalanobis distance in the original space.
    S_inv <- MASS::ginv(stats::cov(reference))
    eig <- eigen(S_inv, symmetric = TRUE)
    W <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)
    reference = reference %*% W
    if (!is.null(query)) {
      query = query %*% W
    }
    dist_fun = "euclidean"
  }


  # calculate the distance matrix
  if (is.null(reference)) {
    dists <- philentropy::distance(query, method = dist_fun)
    if (length(dists) == 1) return(dists)
    diag(dists) <- NA # Exclude self-distance
    } else if (is.null(query)) {
    dists <- philentropy::distance(reference, method = dist_fun)
    if (length(dists) == 1) return(dists)
    diag(dists) <- NA # Exclude self-distance
    } else {
    dists <- philentropy::dist_many_many(query, reference, method = dist_fun)
  }

  if (distance) {
    get_dist <- function(x, k) sort(x)[1:k]
    knn_dists <- t(apply(dists, 1, get_dist, k=k))
    return(knn_dists)
  } else {
    get_index <- function(x, k) order(x)[1:k]
    indices <- t(apply(dists, 1, get_index, k=k))
    return(indices)
  }
}  


.numeric_fct <- function(x) {
  if (is.null(x)) return(NULL)
  catVars <- names(x)[vapply(x, function(z) inherits(z, c("factor", "character")), logical(1))]
  if (length(catVars) == 0) return(x)
  # Convert categorical variables to factor then to numeric
  for (var in catVars) {
    x[[var]] <- as.integer(as.factor(x[[var]]))
  }
  return(x)
}