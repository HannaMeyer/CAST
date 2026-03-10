#' Calculate euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Calculates nearest neighbor distances in geographic space or feature space between training data as well as between prediction locations and training data.
#' Optional, the nearest neighbor distances between test data and training data or between different CV folds is computed.
#' @param x object of class sf, training data locations
#' @param modeldomain SpatRaster, stars or sf object defining the prediction area (see Details)
#' @param dist_space "geographical", "feature" or "time". Should the distance be computed in geographic space, in the normalized multivariate predictor space or in temporal space? (see Details)
#' @param CVtest optional. list or vector. Either a list where each element contains the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm
#' @param CVtrain optional.  A list, where each element contains the data points used for training during the cross validation iteration.
#' Only accepted when CVtest has been specified, and only required if CVtrain is not the opposite of CVtest.
#' Relevant if some data points are excluded, e.g. when using \code{\link{nndm}}.
#' @param cvtrain depreceated.
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
#' @param dist_fun character. Currently covers `euclidean` (default), `gower`, `mahalanobis` and `great_circle`.
#' `gower` and `mahalanobis` only work with `dist_space`="feature", while `great_circle` only works with `dist_space`="geographical". 
#' `mahalanobis` takes into account correlation between predictor values. While `euclidean` and `mahalanobis` only work with numerical variables, 
#' `gower` also works with mixed data including numerical and categorical variables.
#' For `dist_space`="time", currently only the absolute difference is implemented (which equals euclidean distance in one dimension, and is thus covered by `dist_fun=euclidean`)
#' For the geographical space, `great_circle` covers lon/lat coordinates, whereas `euclidean` only works with projected coordinates.
#' @param scale_vars boolean. Should variables be scaled? Only for `dist_space`="feature". 
#' Calculating Gower distances already includes scaling, and manually rescale the data is redundant. 
#' For other distances (Mahalanobis, Euclidean), scaling the data is important. Thus, TRUE by default.
#' @param cvtrain depreceated. Use `CVtrain` instead.
#' @param cvfolds depreceated. Use `CVtest` instead.
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
  cvfolds = NULL
){

  # 1. Input validation & normalization ----------

  ## Guard against old parameter names / depreceated parameters
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


  ## Check that dist_space was correctly defined
  if (dist_space == "geo") dist_space <- "geographical"
  if (!dist_space %in% c("geographical", "feature", "time")) {
    stop("dist_space must be one of 'geographical', 'feature' or 'time'")
  }

  if(dist_space == "time" && dist_fun != "euclidean") stop("Temporal space only supports euclidean distances.")
  if(dist_space == "feature" && dist_fun == "great_circle") stop("Great-circle distances only work with in geographical space.")
  if(dist_space == "geographical" && dist_fun %in% c("mahalanobis", "gower")) stop("Mahalanobis and Gower distances only work in feature space.")

  ## CVtrain
  if(!is.list(CVtest) && !is.null(CVtrain)) stop("CVtrain can only be used when CVtest is supplied as a list.")
  
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
          } else {
            x_cat <- x[,catVars,drop=FALSE]
            pred_points_cat <- pred_points[,catVars,drop=FALSE]
        
            num_cols <- !(names(x) %in% catVars)
            x_num <- x[,num_cols,drop=FALSE]
            pred_points_num <- pred_points[,-which(names(pred_points)%in%catVars),drop=FALSE]

            scale_attr <- attributes(scale(x_num))
            x <- scale(x_num) |> as.data.frame()
            pred_points <- scale(pred_points_num,center=scale_attr$`scaled:center`,
                                scale=scale_attr$`scaled:scale`) |>
              as.data.frame()
        
            x <- as.data.frame(cbind(x, lapply(x_cat, as.factor)))
            pred_points <- as.data.frame(cbind(pred_points, lapply(pred_points_cat, as.factor)))
          
            if (!is.null(testdata)) {
              testdata_cat <- testdata[,catVars,drop=FALSE]
              testdata_num <- testdata[,num_cols]
              testdata <- scale(testdata_num,center=scale_attr$`scaled:center`,
                                scale=scale_attr$`scaled:scale`) |>
               as.data.frame()
        
              testdata <- as.data.frame(cbind(testdata, lapply(testdata_cat, as.factor)))
            }

      }
    }
    
  }

  # 5. Distance computation ----------

  # check that distance function was correctly specified:
  if(!dist_fun %in% c("euclidean", "mahalanobis", "gower", "great_circle")) {
    stop("dist_fun must be one of 'euclidean', 'mahalanobis', 'gower', 'great_circle'")
  }
  if(dist_space == "feature") {
    # Check that dist_fun is appropriate for feature types
    if (!is.null(catVars) && dist_fun != "gower") {
        message("Only gower distances work with categorical features. 'dist_type' was set to 'gower'")
        dist_fun <- "gower"
    } 
  } else if(dist_space == "geographical") {
    if(isTRUE(islonglat) && dist_fun != "great_circle")  {
      message("Only great-circle distances are allowed for lon/lat coordinates. 'dist_fun' was set to 'great_circle'")
      dist_fun <- "great_circle"
    }
  }

  # Calculate NNDs between training points
  s2s <- compute_NND(
    x=x, y=NULL, dist_space = dist_space, dist_type_label = "sample-to-sample",  dist_fun = dist_fun,
    folds = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
  )

  # Calculate NNDs between prediction and training points
  p2s <- compute_NND(
    x=x, y=pred_points, dist_space = dist_space, dist_type_label = "prediction-to-sample",  dist_fun = dist_fun,
    folds = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
  )
  dists <- rbind(s2s, p2s)

  # Calculate NNDs between test points and training points
  if (!is.null(testdata)) {
    t2s <- compute_NND(
      x=x, y=testdata, dist_space = dist_space, dist_type_label = "test-to-sample",  dist_fun = dist_fun,
      folds = NULL, CVtrain = NULL, time_var = time_var, time_unit = time_unit, algorithm = algorithm
    )
    dists <- rbind(dists, t2s)
  }
  # Calculate NNDs between CV folds
  if (!is.null(CVtest)) {
    cvdist <- compute_NND(
      x=x, y=NULL, dist_space = dist_space, dist_type_label = "CV-distances",  dist_fun = dist_fun,
      folds = CVtest, CVtrain = CVtrain, time_var = time_var, time_unit = time_unit, algorithm = algorithm
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
# include the possibility to input train and test folds, not only a vector of fold indices as currently in CV
# OR: allow for NAs in this vector
compute_NND <- function(x, y = NULL, dist_space = c("geographical","feature","time"), 
                      dist_type_label = "sample-to-sample", 
                      dist_fun = c("euclidean","mahalanobis","gower","great_circle"), 
                      folds = NULL, CVtrain = NULL, time_var = NULL, time_unit = "days", algorithm) {
  
  dist_space <- match.arg(dist_space)
  dist_fun <- match.arg(dist_fun)

  # Stop when encountering invalid parameter combinations
  if(!is.null(y) && !is.null(folds)) stop("Currently, `folds` and `y` cannot be specified simultaniously.")

  # Extract coordinates from x (and y)
  coords_x <- if(dist_space == "geographical") sf::st_coordinates(x)[,1:2] else NULL
  coords_y <- if(dist_space == "geographical" && !is.null(y)) sf::st_coordinates(y)[,1:2] else NULL

  # Convert folds to vector
  if(!is.null(folds)){
    if(is.list(folds)) {
      clust <- rep(NA_integer_, nrow(x))
      for(k in seq_along(folds)) clust[folds[[k]]] <- k
    } else clust <- folds
  }

  # Calculate nearest neighbor distances in geographical space
  if(dist_space == "geographical") {
    if(dist_fun == "great_circle") {
      distmat <- if(is.null(y)) sf::st_distance(x) else sf::st_distance(y, x)
      units(distmat) <- NULL
      if(is.null(y)) diag(distmat) <- NA
      min_d <- if(is.null(folds)) {
        apply(distmat, 1, min, na.rm=TRUE)
      } else {
        distclust_distmat(distmat, folds = clust, CVtrain = CVtrain)
      } 
    } else {
      min_d <- if(is.null(y)) {
        if(is.null(folds)) {
          FNN::knn.dist(coords_x, k = 1, algorithm=algorithm)
        } else {
          distclust_euclidean(coords_x, folds = clust, CVtrain = CVtrain, algorithm=algorithm)
        }
      } else {
        FNN::knnx.dist(query = coords_y, data = coords_x, k = 1, algorithm=algorithm)
      }
    }
    
  } else if(dist_space == "feature"){
    # Calculate nearest neighbor distances in feature space
    df_x <- sf::st_drop_geometry(x)
    df_y <- if(!is.null(y)) sf::st_drop_geometry(y) else df_x
    mat_x <- as.matrix(df_x)
    mat_y <- as.matrix(df_y)
    
    if(dist_fun == "mahalanobis"){
      
      if(is.null(folds)) {
  
        S <- if(ncol(mat_x) == 1) matrix(stats::var(mat_x),1,1) else stats::cov(mat_x)
        S_inv <- MASS::ginv(S)
        
        min_d <- vapply(seq_len(nrow(mat_y)), function(i){
          
          idx <- if(is.null(y)) seq_len(nrow(mat_x))[-i] else seq_len(nrow(mat_x))
          
          min(vapply(idx, function(j){
            diff <- mat_y[i,] - mat_x[j,]
            as.numeric(sqrt(t(diff) %*% S_inv %*% diff))
          }, numeric(1)))
          
        }, numeric(1))
        
      } else {
        min_d <- distclust_MD(mat_x, folds = clust, CVtrain = CVtrain)
      }
      
    } else if(dist_fun == "euclidean"){
      if(is.null(folds)) {
        min_d <- if(is.null(y)) FNN::knn.dist(mat_x, k=1, algorithm=algorithm) else
               FNN::knnx.dist(query = mat_y, data = mat_x, k=1, algorithm=algorithm)
      } else {
        min_d <- distclust_euclidean(mat_x, folds = clust, CVtrain = CVtrain, algorithm=algorithm)
      }
      
    } else if(dist_fun == "gower"){
      if(is.null(folds)) {
        min_d <- if(is.null(y)) {
          vapply(1:nrow(df_x), function(i) gower::gower_topn(df_x[i,,drop=FALSE], df_x[-i,,drop=FALSE], n=1)$distance[[1]], numeric(1))
        } else {
          gower_res <- gower::gower_topn(df_y, df_x, n = 1)
          vapply(gower_res$distance, function(x) x[1], numeric(1))
        }
      } else {
        min_d <- distclust_gower(df_x, folds = clust, CVtrain = CVtrain)
      }      
    }
    
  } else if(dist_space == "time"){
    # Calculate nearest neighbor distances in temporal space
    if(dist_fun == "euclidean") {
      time_x <- sf::st_drop_geometry(x)[,time_var]
      if(is.null(folds)) {
        time_y <- if(!is.null(y)) sf::st_drop_geometry(y)[,time_var] else time_x
        dmat <- abs(outer(time_y, time_x, FUN = function(a,b) as.numeric(difftime(a,b,units=time_unit))))
        if(is.null(y)) diag(dmat) <- Inf
        min_d <- apply(dmat, 1, min)
      } else {
        min_d <- distclust_time(time_x, clust, time_unit = time_unit, CVtrain = CVtrain)
      }
      
    }
    
  }
  
  data.frame(dist = min_d,
             what = factor(dist_type_label),
             dist_type = dist_space)
}


## CV-distance functions ----------
# Helper function: Compute out-of-fold NN distance based on a distance matrix (geographical coordinates / numerical variables)
distclust_distmat <- function(distm, folds, CVtrain = NULL){
  alldist <- rep(NA, length(folds))
  
  for(f in unique(folds)){

    test_idx <- which(folds == f)

    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(folds != f)
    }
    
    alldist[test_idx] <- apply(distm[test_idx, train_idx, drop=FALSE], 1, min)
  }
  
  alldist
}


# Helper function: Compute out-of-fold NN distance using euclidean distances and FNN (projected coordinates / numerical variables)
distclust_euclidean <- function(tr_coords, folds, CVtrain = NULL, algorithm){
  alldist <- rep(NA, length(folds))
  for(f in unique(folds)){
    
    test_idx <- which(folds == f)
    
    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(folds != f)
    }

    tr_train <- tr_coords[train_idx,, drop=FALSE]
    tr_test <- tr_coords[test_idx,, drop=FALSE]
    
    alldist[test_idx] <- FNN::knnx.dist(query = tr_test,
                                  data = tr_train,
                                  k = 1,
                                  algorithm = algorithm)
  }
  alldist
}

# Helper function: Compute out-of-fold NN distance using Gower distances (categorical (and numerical) variables)
distclust_gower <- function(tr_coords, folds, CVtrain = NULL){
  alldist <- rep(NA, length(folds))
  
  for(f in unique(folds)){
    
    test_idx <- which(folds == f)
    
    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(folds != f)
    }

    alldist[test_idx] <- c(gower::gower_topn(tr_coords[test_idx,,drop=FALSE],
                                   tr_coords[train_idx,,drop=FALSE],
                                   n = 1))$distance[[1]]
  }
  
  alldist
}

# Helper function: Compute out-of-fold NN distance using Mahalanobis distances (numerical variables)
distclust_MD <- function(tr_coords, folds, CVtrain = NULL){
  tr_mat <- as.matrix(tr_coords)
  S_inv <- MASS::ginv(stats::cov(tr_mat))
  
  alldist <- rep(NA, length(folds))
  
  for(f in unique(folds)){
    test_idx <- which(folds == f)
    
    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(folds != f)
    }
    
    alldist[test_idx] <- apply(tr_mat[test_idx,,drop=FALSE], 1, function(y){
      min(apply(tr_mat[train_idx,,drop=FALSE], 1, function(x){
        sqrt(t(y - x) %*% S_inv %*% (y - x))
      }))
    })
  }
  
  alldist
}

# Helper function: Compute out-of-fold NN distance using absolute time difference (same as euclidean distances in 1D)
# (Only used in geodist at the moment, might be used in future implementations of kNNDM in the temporal space)
distclust_time <- function(time_values, folds, CVtrain = NULL, time_unit = "days"){
  alldist <- rep(NA, length(folds))
  
  for(f in unique(folds)){
    test_idx <- which(folds == f)
    
    if(!is.null(CVtrain)) {
      train_idx <- CVtrain[[f]]
    } else {
      train_idx <- which(folds != f)
    }
    
    diffs <- outer(time_values[test_idx], time_values[train_idx], function(x,y) abs(as.numeric(difftime(x, y, units=time_unit))))
    alldist[test_idx] <- apply(diffs, 1, min)
  }
  
  alldist
}


# Helper function: sample prediction points from the prediction area ----------
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




