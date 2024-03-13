#' Calculate euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Calculates nearest neighbor distances in geographic space or feature space between training data as well as between training data and prediction locations.
#' Optional, the nearest neighbor distances between training data and test data or between training data and CV iterations is computed.
#' @param x object of class sf, training data locations
#' @param modeldomain SpatRaster, stars or sf object defining the prediction area (see Details)
#' @param type "geo" or "feature". Should the distance be computed in geographic space or in the normalized multivariate predictor space (see Details)
#' @param cvfolds optional. list or vector. Either a list where each element contains the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm
#' @param cvtrain optional. List of row indices of x to fit the model to in each CV iteration. If cvtrain is null but cvfolds is not, all samples but those included in cvfolds are used as training data
#' @param testdata optional. object of class sf: Point data used for independent validation
#' @param preddata optional. object of class sf: Point data indicating the locations within the modeldomain to be used as target prediction points. Useful when the prediction objective is a subset of
#' locations within the modeldomain rather than the whole area.
#' @param samplesize numeric. How many prediction samples should be used?
#' @param sampling character. How to draw prediction samples? See \link[sp]{spsample}. Use sampling = "Fibonacci" for global applications.
#' @param variables character vector defining the predictor variables used if type="feature. If not provided all variables included in modeldomain are used.
#' @param timevar optional. character. Column that indicates the date. Only used if type="time".
#' @param time_unit optional. Character. Unit for temporal distances See ?difftime.Only used if type="time".
#' @return A data.frame containing the distances. Unit of returned geographic distances is meters. attributes contain W statistic between prediction area and either sample data, CV folds or test data. See details.
#' @details The modeldomain is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by samplesize) from the spatial extent.
#'     If type = "feature", the argument modeldomain (and if provided then also the testdata and/or preddata) has to include predictors. Predictor values for x, testdata and preddata are optional if modeldomain is a raster.
#'     If not provided they are extracted from the modeldomain rasterStack. If some predictors are categorical (i.e., of class factor or character), gower distances will be used.
#'     W statistic describes the match between the distributions. See Linnenbrink et al (2023) for further details.
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
#'                 type = "feature",
#'                 variables = c("bio_1","bio_12", "elev"))
#' plot(dist)
#'
#' dist <- geodist(x = splotdata[splotdata$Country != "Chile",],
#'                 modeldomain = predictors, cvfolds = folds,
#'                 testdata = splotdata[splotdata$Country == "Chile",],
#'                 type = "feature",
#'                 variables=c("bio_1","bio_12", "elev"))
#' plot(dist)
#'
#'############Distances in temporal space
#' library(lubridate)
#' library(ggplot2)
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' dat <- st_as_sf(dat,coords=c("Easting","Northing"))
#' st_crs(dat) <- 26911
#' trainDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2010,]
#' predictionDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2011,]
#' trainDat$week <- lubridate::week(trainDat$Date)
#' cvfolds <- CreateSpacetimeFolds(trainDat,timevar = "week")
#'
#' dist <- geodist(trainDat,preddata = predictionDat,cvfolds = cvfolds$indexOut,type="time",time_unit="days")
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

geodist <- function(x,
                    modeldomain=NULL,
                    type = "geo",
                    cvfolds=NULL,
                    cvtrain=NULL,
                    testdata=NULL,
                    preddata=NULL,
                    samplesize=2000,
                    sampling = "regular",
                    variables=NULL,
                    timevar=NULL,
                    time_unit="auto"){

  # input formatting ------------
  if(is.null(modeldomain)&!is.null(preddata)){
    modeldomain <- sf::st_bbox(preddata)
  }
  if (inherits(modeldomain, "Raster")) {
    modeldomain <- methods::as(modeldomain,"SpatRaster")
  }
  if (inherits(modeldomain, "stars")) {
    if (!requireNamespace("stars", quietly = TRUE))
      stop("package stars required: install that first")
    modeldomain <- methods::as(modeldomain, "SpatRaster")
  }




  x <- sf::st_transform(x,4326)
  if(type == "feature"){
    if(is.null(variables)){
      variables <- names(modeldomain)
    }
    if(any(!variables%in%names(x))){ # extract variable values of raster:
      message("features are extracted from the modeldomain")
      x <- sf::st_transform(x,sf::st_crs(modeldomain))

      if(class(x)[1]=="sfc_POINT"){
        x <- sf::st_as_sf(x)
      }
      x <- sf::st_as_sf(terra::extract(modeldomain, x, na.rm=FALSE,bind=TRUE))
      x <- sf::st_transform(x,4326)
    }
    if(!is.null(testdata)){
      if(any(!variables%in%names(testdata))){# extract variable values of raster:
        testdata <- sf::st_transform(testdata,sf::st_crs(modeldomain))
        testdata <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(testdata), na.rm=FALSE,bind=TRUE))

        if(any(is.na(testdata))){
          testdata <- na.omit(testdata)
          message("some test data were removed because of NA in extracted predictor values")
        }

        testdata <- sf::st_transform(testdata,4326)
      }
    }
    if(!is.null(preddata)){
      if(any(!variables%in%names(preddata))){# extract variable values of raster:
        preddata <- sf::st_transform(preddata,sf::st_crs(modeldomain))
        preddata <- sf::st_as_sf(terra::extract(modeldomain, terra::vect(preddata), na.rm=FALSE,bind=TRUE))

        if(any(is.na(preddata))){
          preddata <- na.omit(preddata)
          message("some prediction data were removed because of NA in extracted predictor values")
        }

        preddata <- sf::st_transform(preddata,4326)
      }
    }
    # get names of categorical variables
    catVars <- names(x[,variables])[which(sapply(x[,variables], class)%in%c("factor","character"))]
    if(length(catVars)==0) {
      catVars <- NULL
    }
    if(!is.null(catVars)) {
      message(paste0("variable(s) '", catVars, "' is (are) treated as categorical variables"))
    }
  }
  if(type != "feature") {
    catVars <- NULL
  }
  if (type=="time" & is.null(timevar)){
    timevar <- names(which(sapply(x, lubridate::is.Date)))
    message("time variable that has been selected: ",timevar)
  }
  if (type=="time"&time_unit=="auto"){
    time_unit <- units(difftime(sf::st_drop_geometry(x)[,timevar],
                                sf::st_drop_geometry(x)[,timevar]))
  }



  # required steps ----

  ## Sample prediction location from the study area if preddata not available:
  if(is.null(preddata)){
    modeldomain <- sampleFromArea(modeldomain, samplesize, type,variables,sampling, catVars)
    } else{
    modeldomain <- preddata
  }

  # always do sample-to-sample and sample-to-prediction
  s2s <- sample2sample(x, type,variables,time_unit,timevar, catVars)
  s2p <- sample2prediction(x, modeldomain, type, samplesize,variables,time_unit,timevar, catVars)

  dists <- rbind(s2s, s2p)

  # optional steps ----
  ##### Distance to test data:
  if(!is.null(testdata)){
    s2t <- sample2test(x, testdata, type,variables,time_unit,timevar, catVars)
    dists <- rbind(dists, s2t)
  }

  ##### Distance to CV data:
  if(!is.null(cvfolds)){

    cvd <- cvdistance(x, cvfolds, cvtrain, type, variables,time_unit,timevar, catVars)
    dists <- rbind(dists, cvd)
  }
  class(dists) <- c("geodist", class(dists))
  attr(dists, "type") <- type

  if(type=="time"){
    attr(dists, "unit") <- time_unit
  }


  ##### Compute W statistics
  W_sample <- twosamples::wass_stat(dists[dists$what == "sample-to-sample", "dist"],
                                    dists[dists$what == "prediction-to-sample", "dist"])
  attr(dists, "W_sample") <- W_sample
  if(!is.null(testdata)){
    W_test <- twosamples::wass_stat(dists[dists$what == "test-to-sample", "dist"],
                                    dists[dists$what == "prediction-to-sample", "dist"])
    attr(dists, "W_test") <- W_test
  }
  if(!is.null(cvfolds)){
    W_CV <- twosamples::wass_stat(dists[dists$what == "CV-distances", "dist"],
                                  dists[dists$what == "prediction-to-sample", "dist"])
    attr(dists, "W_CV") <- W_CV
  }

  return(dists)
}




# Sample to Sample Distance

sample2sample <- function(x, type,variables,time_unit,timevar, catVars){
  if(type == "geo"){
    sf::sf_use_s2(TRUE)
    d <- sf::st_distance(x)
    diag(d) <- Inf
    min_d <- apply(d, 1, min)
    sampletosample <- data.frame(dist = min_d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "geo")
  }else if(type == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)

    if(!is.null(catVars)) {
      x_cat <- x[,catVars,drop=FALSE]
      x_num <- x[,-which(names(x)%in%catVars),drop=FALSE]
      scaleparam <- attributes(scale(x_num))
      x_num <- data.frame(scale(x_num))
      x <- as.data.frame(cbind(x_num, lapply(x_cat, as.factor)))
      x_clean <- x[complete.cases(x),]
    } else {
      scaleparam <- attributes(scale(x))
      x <- data.frame(scale(x))
      x_clean <- data.frame(x[complete.cases(x),])
    }

    # sample to sample feature distance
    d <- c()
    for (i in 1:nrow(x_clean)){

      if(is.null(catVars)) {
        trainDist <-  FNN::knnx.dist(x_clean[i,],x_clean,k=1)
      } else {
        trainDist <- gower::gower_dist(x_clean[i,],x_clean)
      }

      trainDist[i] <- NA
      d <- c(d,min(trainDist,na.rm=T))
    }
    sampletosample <- data.frame(dist = d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "feature")

  }else if(type == "time"){ # calculate temporal distance matrix
    d <- matrix(ncol=nrow(x),nrow=nrow(x))
    for (i in 1:nrow(x)){
      d[i,] <- abs(difftime(sf::st_drop_geometry(x)[,timevar],
                            sf::st_drop_geometry(x)[i,timevar],
                            units=time_unit))
    }
    diag(d) <- Inf
    min_d <- apply(d, 1, min)
    sampletosample <- data.frame(dist = min_d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "time")
  }
  return(sampletosample)
}


# Sample to Prediction
sample2prediction = function(x, modeldomain, type, samplesize,variables,time_unit,timevar, catVars){

  if(type == "geo"){
    modeldomain <- sf::st_transform(modeldomain, sf::st_crs(x))
    sf::sf_use_s2(TRUE)
    d0 <- sf::st_distance(modeldomain, x)
    min_d0 <- apply(d0, 1, min)
    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("prediction-to-sample"),
                                     dist_type = "geo")

  }else if(type == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    modeldomain <- modeldomain[,variables]
    modeldomain <- sf::st_drop_geometry(modeldomain)

    if(!is.null(catVars)) {

      x_cat <- x[,catVars,drop=FALSE]
      x_num <- x[,-which(names(x)%in%catVars),drop=FALSE]
      scaleparam <- attributes(scale(x_num))
      x_num <- data.frame(scale(x_num))

      modeldomain_num <- modeldomain[,-which(names(modeldomain)%in%catVars),drop=FALSE]
      modeldomain_cat <- modeldomain[,catVars,drop=FALSE]
      modeldomain_num <- data.frame(scale(modeldomain_num,center=scaleparam$`scaled:center`,
                                      scale=scaleparam$`scaled:scale`))

      x <- as.data.frame(cbind(x_num, lapply(x_cat, as.factor)))
      x_clean <- x[complete.cases(x),]
      modeldomain <- as.data.frame(cbind(modeldomain_num, lapply(modeldomain_cat, as.factor)))

    } else {
      scaleparam <- attributes(scale(x))
      x <- data.frame(scale(x))
      x_clean <- x[complete.cases(x),]

      modeldomain <- data.frame(scale(modeldomain,center=scaleparam$`scaled:center`,
                                      scale=scaleparam$`scaled:scale`))
    }


    target_dist_feature <- c()
    for (i in 1:nrow(modeldomain)){

      if(is.null(catVars)) {
        trainDist <-  FNN::knnx.dist(modeldomain[i,],x_clean,k=1)
      } else {
        trainDist <- gower::gower_dist(modeldomain[i,], x_clean)
      }

      target_dist_feature <- c(target_dist_feature,min(trainDist,na.rm=T))
    }
    sampletoprediction <- data.frame(dist = target_dist_feature,
                                     what = "prediction-to-sample",
                                     dist_type = "feature")
  }else if(type == "time"){

    min_d0 <- c()
    for (i in 1:nrow(modeldomain)){
      min_d0[i] <- min(abs(difftime(sf::st_drop_geometry(modeldomain)[i,timevar],
                                    sf::st_drop_geometry(x)[,timevar],
                                    units=time_unit)))
    }

    sampletoprediction <- data.frame(dist = min_d0,
                                     what = factor("prediction-to-sample"),
                                     dist_type = "time")

  }

  return(sampletoprediction)
}


# sample to test


sample2test <- function(x, testdata, type,variables,time_unit,timevar, catVars){

  if(type == "geo"){
    testdata <- sf::st_transform(testdata,4326)
    d_test <- sf::st_distance(testdata, x)
    min_d_test <- apply(d_test, 1, min)

    dists_test <- data.frame(dist = min_d_test,
                             what = factor("test-to-sample"),
                             dist_type = "geo")


  }else if(type == "feature"){


    x <- x[,variables]
    x <- sf::st_drop_geometry(x)
    testdata <- testdata[,variables]
    testdata <- sf::st_drop_geometry(testdata)

    if(!is.null(catVars)) {

      x_cat <- x[,catVars,drop=FALSE]
      x_num <- x[,-which(names(x)%in%catVars),drop=FALSE]
      scaleparam <- attributes(scale(x_num))
      x_num <- data.frame(scale(x_num))

      testdata_num <- testdata[,-which(names(testdata)%in%catVars),drop=FALSE]
      testdata_cat <- testdata[,catVars,drop=FALSE]
      testdata_num <- data.frame(scale(testdata_num,center=scaleparam$`scaled:center`,
                                          scale=scaleparam$`scaled:scale`))

      x <- as.data.frame(cbind(x_num, lapply(x_cat, as.factor)))
      x_clean <- x[complete.cases(x),]
      testdata <- as.data.frame(cbind(testdata_num, lapply(testdata_cat, as.factor)))

    } else {
      scaleparam <- attributes(scale(x))
      x <- data.frame(scale(x))
      x_clean <- x[complete.cases(x),]

      testdata <- data.frame(scale(testdata,center=scaleparam$`scaled:center`,
                                      scale=scaleparam$`scaled:scale`))
    }


    test_dist_feature <- c()
    for (i in 1:nrow(testdata)){

      if(is.null(catVars)) {
        testDist <- FNN::knnx.dist(testdata[i,],x_clean,k=1)
      } else {
        testDist <- gower::gower_dist(testdata[i,], x_clean)
      }
      test_dist_feature <- c(test_dist_feature,min(testDist,na.rm=T))
    }
    dists_test <- data.frame(dist = test_dist_feature,
                             what = "test-to-sample",
                             dist_type = "feature")
  }else if (type=="time"){
    min_d0 <- c()
    for (i in 1:nrow(testdata)){
      min_d0[i] <- min(abs(difftime(sf::st_drop_geometry(testdata)[i,timevar],
                                    sf::st_drop_geometry(x)[,timevar],
                                    units=time_unit)))
    }

    dists_test <- data.frame(dist = min_d0,
                             what = factor("test-to-sample"),
                             dist_type = "time")



  }
  return(dists_test)
}



# between folds

cvdistance <- function(x, cvfolds, cvtrain, type, variables,time_unit,timevar, catVars){

  if(!is.null(cvfolds)&!is.list(cvfolds)){ # restructure input if CVtest only contains the fold ID
    tmp <- list()
    for (i in unique(cvfolds)){
      tmp[[i]] <- which(cvfolds==i)
    }
    cvfolds <- tmp
  }


  if(type == "geo"){
    d_cv <- c()
    for (i in 1:length(cvfolds)){

      if(!is.null(cvtrain)){
        d_cv_tmp <- sf::st_distance(x[cvfolds[[i]],], x[cvtrain[[i]],])
      }else{
        d_cv_tmp <- sf::st_distance(x[cvfolds[[i]],], x[-cvfolds[[i]],])
      }
      d_cv <- c(d_cv,apply(d_cv_tmp, 1, min))
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "geo")


  }else if(type == "feature"){
    x <- x[,variables]
    x <- sf::st_drop_geometry(x)

    if(is.null(catVars)) {
      x <- data.frame(scale(x))
    } else {
      x_cat <- x[,catVars,drop=FALSE]
      x_num <- x[,-which(names(x)%in%catVars),drop=FALSE]
      scaleparam <- attributes(scale(x_num))
      x_num <- data.frame(scale(x_num))
      x <- as.data.frame(cbind(x_num, lapply(x_cat, as.factor)))
    }

    d_cv <- c()
    for(i in 1:length(cvfolds)){

      if(!is.null(cvtrain)){
        testdata_i <- x[cvfolds[[i]],]
        traindata_i <- x[cvtrain[[i]],]
      }else{
        testdata_i <- x[cvfolds[[i]],]
        traindata_i <- x[-cvfolds[[i]],]
      }

      testdata_i <- testdata_i[complete.cases(testdata_i),]
      traindata_i <- traindata_i[complete.cases(traindata_i),]

      for (k in 1:nrow(testdata_i)){

        if(is.null(catVars)) {
          trainDist <-  tryCatch(FNN::knnx.dist(testdata_i[k,],traindata_i,k=1),
                                 error = function(e)e)
          if(inherits(trainDist, "error")){
            trainDist <- NA
            message("warning: no distance could be calculated for a fold.
                  Possibly because predictor values are NA")
          }
        } else {
          trainDist <-  tryCatch(gower::gower_dist(testdata_i[i,], traindata_i),
                                 error = function(e)e)
          if(inherits(trainDist, "error")){
            trainDist <- NA
            message("warning: no distance could be calculated for a fold.
                  Possibly because predictor values are NA")
          }
        }


        trainDist[k] <- NA
        d_cv <- c(d_cv,min(trainDist,na.rm=T))
      }
    }

    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "feature")

  }else if(type == "time"){
    d_cv <- c()
    d_cv_tmp <- c()
    for (i in 1:length(cvfolds)){
      if(!is.null(cvtrain)){
        for (k in 1:length(cvfolds[[i]])){
          d_cv_tmp[k] <- min(abs(difftime(sf::st_drop_geometry(x)[cvfolds[[i]][k],timevar],
                                          sf::st_drop_geometry(x)[cvtrain[[i]],timevar],
                                          units=time_unit)))
        }
      }else{
        for (k in 1:length(cvfolds[[i]])){
          d_cv_tmp[k] <- min(abs(difftime(sf::st_drop_geometry(x)[cvfolds[[i]][k],timevar],
                                          sf::st_drop_geometry(x)[-cvfolds[[i]],timevar],
                                          units=time_unit)))
        }
      }
      d_cv <- c(d_cv,d_cv_tmp)
    }


    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "time")

  }

  return(dists_cv)
}





sampleFromArea <- function(modeldomain, samplesize, type,variables,sampling, catVars){

  ##### Distance to prediction locations:
  # regularly spread points (prediction locations):
  # see https://edzer.github.io/OGH21/
  if(inherits(modeldomain, "Raster")){
    modeldomain <- terra::rast(modeldomain)
  }

  if(inherits(modeldomain, "SpatRaster")) {
    if(samplesize>terra::ncell(modeldomain)){
      samplesize <- terra::ncell(modeldomain)
      message(paste0("samplesize for new data shouldn't be larger than number of pixels.
              Samplesize was reduced to ",terra::ncell(modeldomain)))
    }
    #create mask to sample from:
    template <- modeldomain[[1]]
    terra::values(template)[!is.na(terra::values(template))] <-1
    modeldomainextent <- terra::as.polygons(template) |>
      sf::st_as_sf() |>
      sf::st_geometry()
  }else{
    modeldomainextent <- modeldomain
  }

  sf::sf_use_s2(FALSE)
  sf::st_as_sf(modeldomainextent) |>
    sf::st_transform(4326) -> bb

  methods::as(bb, "Spatial") |>
    sp::spsample(n = samplesize, type = sampling)  |>
    sf::st_as_sfc() |>
    sf::st_set_crs(4326) -> predictionloc

  predictionloc <- sf::st_as_sf(predictionloc)


  if(type == "feature"){

    if(is.null(catVars)) {
      modeldomain <- terra::project(modeldomain, "epsg:4326")
    } else {
      modeldomain <- terra::project(modeldomain, "epsg:4326", method="near")
    }

    predictionloc <- sf::st_as_sf(terra::extract(modeldomain,terra::vect(predictionloc),bind=TRUE))
    predictionloc <- na.omit(predictionloc)
  }

  return(predictionloc)

}




