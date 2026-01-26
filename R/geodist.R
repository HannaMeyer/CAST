#' Calculate euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Calculates nearest neighbor distances in geographic space or feature space between training data as well as between prediction locations and training data.
#' Optional, the nearest neighbor distances between test data and training data or between different CV folds is computed.
#' @param x object of class sf, training data locations
#' @param modeldomain SpatRaster, stars or sf object defining the prediction area (see Details)
#' @param space "geographical", "feature" or "time". Should the distance be computed in geographic space, in the normalized multivariate predictor space or in temporal space? (see Details)
#' @param cvfolds optional. list or vector. Either a list where each element contains the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm
#' @param testdata optional. object of class sf: Point data used for independent validation. May already include the predictor values if `space`=feature.
#' @param preddata optional. object of class sf: Point data indicating the locations within the modeldomain to be used as target prediction points. Useful when the prediction objective is a subset of
#' locations within the modeldomain rather than the whole area. May already include the predictor values if `space`="feature".
#' @param samplesize numeric. How many prediction samples should be used?
#' @param sampling character. How to draw prediction samples? See \link[sf]{st_sample} for modeldomains that are sf objects and \link[terra]{spatSample} for raster objects.
#' Use sampling = "Fibonacci" for global applications (raster objects will be transformed to polygons in this case).
#' @param variables character vector defining the predictor variables used if space="feature". If not provided all variables included in modeldomain are used.
#' @param timevar optional. character. Column that indicates the date. Only used if space="time".
#' @param time_unit optional. Character. Unit for temporal distances See ?difftime.Only used if space="time".
#' @param algorithm see \code{\link[FNN]{knnx.dist}} and \code{\link[FNN]{knnx.index}}
#' @param useMD boolean. Only for `space`=feature: shall the Mahalanobis distance be calculated instead of Euclidean? Only works with numerical variables.
#' @return A data.frame containing the distances. Unit of returned geographic distances is meters. attributes contain W statistic between prediction area and either sample data, CV folds or test data. See details.
#' @details The modeldomain is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by samplesize) from the spatial extent (if no `preddata` are supplied).
#'     If `space` = "feature", the argument modeldomain has to be a raster and include predictors. The only exception is when the provided training data and preddata already include the predictor values.
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
#' dist <- geodist(x=splotdata, modeldomain=studyArea, cvfolds=folds)
#' # Using density functions
#' plot(dist)
#' # Using ECDFs (relevant for nndm and knnmd methods)
#' plot(dist, stat="ecdf")
#'
#' ########### Distances in the feature space:
#' predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))
#' dist <- geodist(x = splotdata,
#'                 modeldomain = predictors,
#'                 space = "feature",
#'                 variables = c("bio_1","bio_12", "elev"))
#' plot(dist)
#'
#' dist <- geodist(x = splotdata[splotdata$Country != "Chile",],
#'                 modeldomain = predictors,
#'                 testdata = splotdata[splotdata$Country == "Chile",],
#'                 space = "feature",
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
#' cvfolds <- CreateSpacetimeFolds(trainDat,timevar = "week")
#'
#' dist <- geodist(trainDat,preddata = predictionDat,cvfolds = cvfolds$indexOut,
#'    space="time",time_unit="days")
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
  space = "geographical",
  cvfolds = NULL,
  testdata = NULL,
  preddata = NULL,
  samplesize = 2000,
  sampling = "regular",
  variables = NULL,
  timevar = NULL,
  time_unit = "auto",
  algorithm = "brute",
  useMD = FALSE
){

  # 1. Input validation & normalization ----------

  # Check that space was correctly defined
  if (space == "geo") space <- "geographical"
  if (!space %in% c("geographical", "feature", "time")) {
    stop("Space must be one of 'geographical', 'feature' or 'time'")
  }

  # Check for different raster formats
  if (inherits(modeldomain, "Raster")) {
    modeldomain <- methods::as(modeldomain,"SpatRaster")
  }
  if (inherits(modeldomain, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    modeldomain <- methods::as(modeldomain, "SpatRaster")
  }

  # Check for time variable and unit
  if (space == "time") {
    if (is.null(timevar)) {
      timevar <- names(which(sapply(x, lubridate::is.Date)))
      message("time variable selected: ", timevar)
    }
    if (time_unit == "auto") {
      time_unit <- units(difftime(
        sf::st_drop_geometry(x)[, timevar],
        sf::st_drop_geometry(x)[, timevar]
      ))
    }
  }


  # 2. CRS harmonization ----------

  # Retrieve the CRS of preddata (or modeldomain) and use it as a reference
  ref_crs <- if (!is.null(preddata)) {
    sf::st_crs(preddata)
  } else if (!is.null(modeldomain) && inherits(modeldomain, "sf")) {
    sf::st_crs(modeldomain)
  } else {
    sf::st_crs(x)
  }

  # transform training points to the CRS of preddata (or modeldomain)
  if (!is.na(ref_crs) && !is.na(sf::st_crs(x)) && sf::st_crs(x) != ref_crs) {
    message("Transforming training data CRS")
    x <- sf::st_transform(x, ref_crs)
  }

  # transform test points to the CRS of preddata (or modeldomain)
  if (!is.null(testdata) &&
      !is.na(sf::st_crs(testdata)) &&
      !is.na(sf::st_crs(x)) &&
      sf::st_crs(testdata) != sf::st_crs(x)) {
    message("Transforming test data CRS")
    testdata <- sf::st_transform(testdata, sf::st_crs(x))
  }

  # Check if coordinates are longitutude/latitude
  islonglat <- if (is.na(ref_crs)) {
    warning("Missing CRS of the modeldomain or prediction points. Assuming projected CRS.")
    FALSE
  } else {
    sf::st_is_longlat(ref_crs)
  }

  # extract the coordinates of the training points
  tcoords <- sf::st_coordinates(x)[, 1:2]


  # 3. Prediction point generation ----------

  # Sample prediction points from the study area (only if no preddata are supplied)
  pred_points <- if (is.null(preddata)) {
    sampleFromArea(modeldomain, samplesize, space, variables, sampling)
  } else {
    preddata
  }


  # 4. Feature-space preparation (if needed) ----------

  catVars <- NULL

  if (space == "feature") {

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
      x <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(x), bind = TRUE))
    }
    x <- x[, variables]

    if (!is.null(testdata) && any(!variables %in% names(testdata))) {
      testdata <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(testdata), bind = TRUE))
    }
    if (!is.null(testdata)) testdata <- testdata[, variables]

    if (any(!variables %in% names(pred_points))) {
      pred_points <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(pred_points), bind = TRUE))
    }
    pred_points <- pred_points[, variables]

    # Detect categorical variables
    catVars <- names(x)[sapply(x, function(z) inherits(z, c("factor", "character")))]
    if (length(catVars) == 0) catVars <- NULL

    # Drop the geometry of the training and prediction (and test) points
    x <- sf::st_drop_geometry(x)
    pred_points <- sf::st_drop_geometry(pred_points)
    if (!is.null(testdata)) testdata <- sf::st_drop_geometry(testdata)

    # Scale the training, prediction and test points
    if (is.null(catVars)) {
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

  # Calculate NNDs between training points
  s2s <- sample2sample(
    x, space, variables, time_unit, timevar,
    catVars, algorithm, useMD, islonglat, tcoords
  )

  # Calculate NNDs between prediction and training points
  p2s <- prediction2sample(
    x, pred_points, space, samplesize, variables,
    time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords
  )

  dists <- rbind(s2s, p2s)

  # Calculate NNDs between test points and training points
  if (!is.null(testdata)) {
    dists <- rbind(dists, test2sample(x, testdata, space, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords))
  }

  # Calculate NNDs between CV folds
  if (!is.null(cvfolds)) {
    dists <- rbind(dists, cvdistance(x, cvfolds, space, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords))
  }


  # 6. Post-processing ----------

  class(dists) <- c("geodist", class(dists))
  attr(dists, "space") <- space
  if (space == "time") attr(dists, "unit") <- time_unit

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

  if (!is.null(cvfolds)) {
    attr(dists, "W_CV") <- twosamples::wass_stat(
      dists[dists$what == "CV-distances", "dist"],
      dists[dists$what == "prediction-to-sample", "dist"]
    )
  }

  dists
}



## Function definitions ----------

# Sample to sample distance calculation
sample2sample <- function(x, space, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords){
  if(space == "geographical"){

    if(isTRUE(islonglat)){
      distmat <- sf::st_distance(x)
      units(distmat) <- NULL
      diag(distmat) <- NA
      min_d <- apply(distmat, 1, function(x) min(x, na.rm=TRUE))
    }else{
      min_d <- c(FNN::knn.dist(tcoords, k = 1, algorithm=algorithm))
    }

    sampletosample <- data.frame(dist = min_d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "geographical")
  }else if(space == "feature"){
    
    if(is.null(catVars)) {
      if(isTRUE(useMD)) {
        tpoints_mat <- as.matrix(x)

        # use Mahalanobis distances
        if (dim(tpoints_mat)[2] == 1) {
          S <- matrix(stats::var(tpoints_mat), 1, 1)
          tpoints_mat <- as.matrix(tpoints_mat, ncol = 1)
        } else {
          S <- stats::cov(tpoints_mat)
        }
        S_inv <- MASS::ginv(S)

        # calculate distance matrix
        distmat <- matrix(nrow=nrow(x), ncol=nrow(x))
        distmat <- sapply(1:nrow(distmat), function(i) {
          sapply(1:nrow(distmat), function(j) {
            sqrt(t(tpoints_mat[i,] - tpoints_mat[j,]) %*% S_inv %*% (tpoints_mat[i,] - tpoints_mat[j,]))
          })
        })
        diag(distmat) <- NA

        d <- apply(distmat, 1, min, na.rm=TRUE)
      } else {
        d <- c(FNN::knn.dist(x, k = 1, algorithm=algorithm))
      }
    } else {
      # use Gower distances if categorical variables are present
      d <- sapply(1:nrow(x), function(i) gower::gower_topn(x[i,], x[-i,], n=1)$distance[[1]])
    }

    sampletosample <- data.frame(dist = d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "feature")

  }else if(space == "time"){ # calculate temporal distance matrix
    time_values <- sf::st_drop_geometry(x)[, timevar]
    d <- abs(outer(time_values, time_values, FUN=function(a,b) as.numeric(difftime(a,b,units=time_unit))))
    diag(d) <- Inf
    min_d <- apply(d, 1, min)
    sampletosample <- data.frame(dist = min_d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "time")
  }
  return(sampletosample)
}


# Prediction to sample distance calculation
prediction2sample = function(x, modeldomain, space, samplesize, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords){

  if(space == "geographical"){

    # calculate the NNDs between prediction points and training points
    if(isTRUE(islonglat)){
      d0 <- sf::st_distance(modeldomain, x)
      units(d0) <- NULL
      min_d0 <- apply(d0, 1, min)
    }else{
      min_d0 <- c(FNN::knnx.dist(query = sf::st_coordinates(modeldomain)[,1:2],
                              data = tcoords, k = 1, algorithm=algorithm))
    }

    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("prediction-to-sample"),
                                     dist_type = "geographical")

  }else if(space == "feature"){

    if(is.null(catVars)) {
      
      if(isTRUE(useMD)) {

        tpoints_mat <- as.matrix(x)
        predpoints_mat <- as.matrix(modeldomain)

        # use Mahalanobis distances
        if (dim(tpoints_mat)[2] == 1) {
          S <- matrix(stats::var(tpoints_mat), 1, 1)
          tpoints_mat <- as.matrix(tpoints_mat, ncol = 1)
        } else {
          S <- stats::cov(tpoints_mat)
        }
        S_inv <- MASS::ginv(S)

        target_dist_feature <- sapply(1:dim(predpoints_mat)[1], function(y) {
          min(sapply(1:dim(tpoints_mat)[1], function(x) {
            sqrt(t(predpoints_mat[y,] - tpoints_mat[x,]) %*% S_inv %*% (predpoints_mat[y,] - tpoints_mat[x,]))
          }))
        })
      } else {
        target_dist_feature <- c(FNN::knnx.dist(query = modeldomain, data = x, k = 1, algorithm=algorithm))
      }

    } else {
      target_dist_feature <- c(gower::gower_topn(modeldomain, x, n = 1)$distance)
    }

    sampletoprediction <- data.frame(dist = target_dist_feature,
                                     what = "prediction-to-sample",
                                     dist_type = "feature")
  }else if(space == "time"){

    time_train <- sf::st_drop_geometry(x)[, timevar]
    time_pred <- sf::st_drop_geometry(modeldomain)[, timevar]
    dmat <- abs(outer(time_pred, time_train, FUN=function(a,b) as.numeric(difftime(a,b,units=time_unit))))
    min_d0 <- apply(dmat, 1, min)

    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("prediction-to-sample"),
                                     dist_type = "time")

  }

  return(sampletoprediction)
}


# Test to sample distance calculation
test2sample <- function(x, testdata, space, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords){

  if(space == "geographical"){

    # calculate the NNDs between test points and training points
    if(isTRUE(islonglat)){
      d_test <- sf::st_distance(testdata, x)
      units(d_test) <- NULL
      min_d_test <- apply(d_test, 1, min)
    }else{
      min_d_test <- c(FNN::knnx.dist(query = sf::st_coordinates(testdata)[,1:2],
                              data = tcoords, k = 1, algorithm=algorithm))
    }
  
   dists_test <- data.frame(dist = min_d_test,
                             what = factor("test-to-sample"),
                             dist_type = "geographical")


  }else if(space == "feature"){

    if(is.null(catVars)) {
      
      if(isTRUE(useMD)) {

        tpoints_mat <- as.matrix(x)
        testpoints_mat <- as.matrix(testdata)

        # use Mahalanobis distances
        if (dim(tpoints_mat)[2] == 1) {
          S <- matrix(stats::var(tpoints_mat), 1, 1)
          tpoints_mat <- as.matrix(tpoints_mat, ncol = 1)
        } else {
          S <- stats::cov(tpoints_mat)
        }
        S_inv <- MASS::ginv(S)

        test_dist_feature <- sapply(1:dim(testpoints_mat)[1], function(y) {
          min(sapply(1:dim(tpoints_mat)[1], function(x) {
            sqrt(t(testpoints_mat[y,] - tpoints_mat[x,]) %*% S_inv %*% (testpoints_mat[y,] - tpoints_mat[x,]))
          }))
        })
      } else {
        test_dist_feature <- c(FNN::knnx.dist(query = testdata, data = x, k = 1, algorithm=algorithm))
      }

    } else {
      test_dist_feature <- c(gower::gower_topn(testdata, x, n = 1)$distance)
    }

    dists_test <- data.frame(dist = test_dist_feature,
                             what = "test-to-sample",
                             dist_type = "feature")
  }else if (space=="time"){
    time_train <- sf::st_drop_geometry(x)[, timevar]
    time_test <- sf::st_drop_geometry(testdata)[, timevar]
    dmat <- abs(outer(time_test, time_train, FUN=function(a,b) as.numeric(difftime(a,b,units=time_unit))))
    min_d0 <- apply(dmat, 1, min)

    dists_test <- data.frame(dist = min_d0,
                             what = factor("test-to-sample"),
                             dist_type = "time")



  }
  return(dists_test)
}


# Calculate distances between folds
cvdistance <- function(x, cvfolds, space, variables, time_unit, timevar, catVars, algorithm, useMD, islonglat, tcoords){
  
  # Convert cvfold list to vector
  if(is.list(cvfolds)){
    n <- max(unlist(cvfolds))
    clust <- integer(n)

    for (k in seq_along(cvfolds)) {
      clust[cvfolds[[k]]] <- k
    }
  } else {
    clust <- cvfolds
  }

  if(space == "geographical"){

    if(isTRUE(islonglat)){
      distmat <- sf::st_distance(x)
      units(distmat) <- NULL
      diag(distmat) <- NA
      d_cv <- distclust_distmat(distmat, clust)
    }else{
      d_cv <- distclust_euclidean(tcoords, clust, algorithm=algorithm)
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "geographical")


  }else if(space == "feature"){

    if(is.null(catVars)) {
      if(isTRUE(useMD)) {
        d_cv <- distclust_MD(x, clust)
      } else {
        d_cv <- distclust_euclidean(x, clust, algorithm=algorithm)
      }

    } else {
      d_cv <- distclust_gower(x, clust)
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "feature")

  }else if(space == "time"){
    time_values <- sf::st_drop_geometry(x)[, timevar]
    d_cv <- distclust_time(time_values, clust, time_unit = time_unit)

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "time")

  }

  return(dists_cv)
}


# Sample prediction points from the prediction area
sampleFromArea <- function(modeldomain, samplesize, space, variables, sampling){

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
    predictionloc <- sf::st_sample(template, size = samplesize, type = sampling) |> 
      sf::st_set_crs(sf::st_crs(modeldomain))

  } else {
    # Sample points from a Polygon
    message(paste0("Sampling ", samplesize, " prediction locations from the modeldomain vector."))
    predictionloc <- sf::st_sample(modeldomain, size = samplesize, type = sampling) |> 
      sf::st_set_crs(sf::st_crs(modeldomain))
  }
  return(predictionloc)

}




