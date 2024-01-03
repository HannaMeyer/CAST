#' Plot euclidean nearest neighbor distances in geographic space or feature space
#'
#' @description Density plot of nearest neighbor distances in geographic space or feature space between training data as well as between training data and prediction locations.
#' Optional, the nearest neighbor distances between training data and test data or between training data and CV iterations is shown.
#' The plot can be used to check the suitability of a chosen CV method to be representative to estimate map accuracy. Alternatively distances can also be calculated in the multivariate feature space.
#'
#' @param x object of class sf, training data locations
#' @param modeldomain SpatRaster, stars or sf object defining the prediction area (see Details)
#' @param type "geo" or "feature". Should the distance be computed in geographic space or in the normalized multivariate predictor space (see Details)
#' @param cvfolds optional. list or vector. Either a list where each element contains the data points used for testing during the cross validation iteration (i.e. held back data).
#' Or a vector that contains the ID of the fold for each training point. See e.g. ?createFolds or ?CreateSpacetimeFolds or ?nndm
#' @param cvtrain optional. List of row indices of x to fit the model to in each CV iteration. If cvtrain is null but cvfolds is not, all samples but those included in cvfolds are used as training data
#' @param testdata optional. object of class sf: Data used for independent validation
#' @param samplesize numeric. How many prediction samples should be used?
#' @param sampling character. How to draw prediction samples? See \link[sp]{spsample}. Use sampling = "Fibonacci" for global applications.
#' @param variables character vector defining the predictor variables used if type="feature. If not provided all variables included in modeldomain are used.
#' @param unit character. Only if type=="geo" and only applied to the plot. Supported: "m" or "km".
#' @param stat "density" for density plot or "ecdf" for empirical cumulative distribution function plot.
#' @param showPlot logical
#' @return A list including the plot and the corresponding data.frame containing the distances. Unit of returned geographic distances is meters.
#' @details The modeldomain is a sf polygon or a raster that defines the prediction area. The function takes a regular point sample (amount defined by samplesize) from the spatial extent.
#'     If type = "feature", the argument modeldomain (and if provided then also the testdata) has to include predictors. Predictor values for x are optional if modeldomain is a raster.
#'     If not provided they are extracted from the modeldomain rasterStack.
#' @note See Meyer and Pebesma (2022) for an application of this plotting function
#' @seealso \code{\link{nndm}}
#' @import ggplot2
#' @author Hanna Meyer, Edzer Pebesma, Marvin Ludwig
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(caret)
#'
#' ########### prepare sample data:
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing")],
#'   by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' pts_train <- pts[1:29,]
#' pts_test <- pts[30:42,]
#' studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
#' studyArea <- studyArea[[c("DEM","TWI", "NDRE.M", "NDRE.Sd", "Bt")]]
#'
#' ########### Distance between training data and new data:
#' dist <- plot_geodist(pts_train,studyArea)
#'
#' ########### Distance between training data, new data and test data:
#' #mapview(pts_train,col.regions="blue")+mapview(pts_test,col.regions="red")
#' dist <- plot_geodist(pts_train,studyArea,testdata=pts_test)
#'
#' ########### Distance between training data, new data and CV folds:
#' folds <- createFolds(1:nrow(pts_train),k=3,returnTrain=FALSE)
#' dist <- plot_geodist(x=pts_train, modeldomain=studyArea, cvfolds=folds)
#'
#' ## or use nndm to define folds
#' AOI <- as.polygons(rast(studyArea), values = F) |>
#'   st_as_sf() |>
#'   st_union() |>
#'   st_transform(crs = st_crs(pts_train))
#' nndm_pred <- nndm(pts_train, AOI)
#' dist <- plot_geodist(x=pts_train, modeldomain=studyArea,
#'     cvfolds=nndm_pred$indx_test, cvtrain=nndm_pred$indx_train)
#'
#' ########### Distances in the feature space:
#' plot_geodist(x=pts_train, modeldomain=studyArea,
#'     type = "feature",variables=c("DEM","TWI", "NDRE.M"))
#'
#' dist <- plot_geodist(x=pts_train, modeldomain=studyArea, cvfolds = folds, testdata = pts_test,
#'     type = "feature",variables=c("DEM","TWI", "NDRE.M"))
#'
#'############ Example for a random global dataset
#'############ (refer to figure in Meyer and Pebesma 2022)
#'library(sf)
#'library(rnaturalearth)
#'library(ggplot2)
#'
#'### Define prediction area (here: global):
#'ee <- st_crs("+proj=eqearth")
#'co <- ne_countries(returnclass = "sf")
#'co.ee <- st_transform(co, ee)
#'
#'### Simulate a spatial random sample
#'### (alternatively replace pts_random by a real sampling dataset (see Meyer and Pebesma 2022):
#'sf_use_s2(FALSE)
#'pts_random <- st_sample(co.ee, 2000, exact=FALSE)
#'
#'### See points on the map:
#'ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
#'      geom_sf(data = pts_random, color = "#F8766D",size=0.5, shape=3) +
#'      guides(fill = FALSE, col = FALSE) +
#'      labs(x = NULL, y = NULL)
#'
#'### plot distances:
#'dist <- plot_geodist(pts_random,co.ee,showPlot=FALSE)
#'dist$plot+scale_x_log10(labels=round)
#'}
#' @export

plot_geodist <- function(x,
                         modeldomain,
                         type = "geo",
                         cvfolds=NULL,
                         cvtrain=NULL,
                         testdata=NULL,
                         samplesize=2000,
                         sampling = "regular",
                         variables=NULL,
                         unit="m",
                         stat = "density",
                         showPlot=TRUE){


  message("plot_geodist() is deprecated and will be removed soon. \n Please use generic plot() function on geodist object from geodist().")


  # input formatting ------------
  if (inherits(modeldomain, "Raster")) {
#    if (!requireNamespace("raster", quietly = TRUE))
#      stop("package raster required: install that first")
    message("Raster will soon not longer be supported. Use terra or stars instead")
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

      #x <- sf::st_as_sf(raster::extract(modeldomain, x, df = TRUE, sp = TRUE))
      x <- sf::st_as_sf(terra::extract(modeldomain, x, na.rm=FALSE,bind=TRUE))
      x <- sf::st_transform(x,4326)
    }
    if(!is.null(testdata)){
      if(any(!variables%in%names(testdata))){# extract variable values of raster:
        testdata <- sf::st_transform(testdata,sf::st_crs(modeldomain))
        #testdata <- sf::st_as_sf(raster::extract(modeldomain, testdata, df = TRUE, sp = TRUE))
        testdata <- sf::st_as_sf(terra::extract(modeldomain, testdata, na.rm=FALSE,bind=TRUE))

        if(any(is.na(testdata))){
          testdata <- na.omit(testdata)
          message("some test data were removed because of NA in extracted predictor values")
        }

        testdata <- sf::st_transform(testdata,4326)
      }
    }
  }


  # required steps ----

  ## Sample prediction location from the study area:
  modeldomain <- sampleFromArea(modeldomain, samplesize, type,variables,sampling)

  # always do sample-to-sample and sample-to-prediction
  s2s <- sample2sample(x, type,variables)
  s2p <- sample2prediction(x, modeldomain, type, samplesize,variables)

  dists <- rbind(s2s, s2p)

  # optional steps ----
  ##### Distance to test data:
  if(!is.null(testdata)){
    s2t <- sample2test(x, testdata, type,variables)
    dists <- rbind(dists, s2t)
  }

  ##### Distance to CV data:
  if(!is.null(cvfolds)){
    cvd <- cvdistance(x, cvfolds, cvtrain, type, variables)
    dists <- rbind(dists, cvd)
  }

  # Compile output and plot data ----
  p <- .plot.nnd(dists,type,unit,stat)

  if(showPlot){
    print(p)
  }

  out <- list(p,dists)
  names(out) <- c("plot","distances")
  return(out)
}




# Sample to Sample Distance

sample2sample <- function(x, type,variables){

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
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))
    x_clean <- data.frame(x[complete.cases(x),])
    # sample to sample feature distance
    d <- c()
    for (i in 1:nrow(x_clean)){

      trainDist <-  FNN::knnx.dist(x_clean[i,],x_clean,k=1)

      trainDist[i] <- NA
      d <- c(d,min(trainDist,na.rm=T))
    }
    sampletosample <- data.frame(dist = d,
                                 what = factor("sample-to-sample"),
                                 dist_type = "feature")

  }
  return(sampletosample)
}


# Sample to Prediction
sample2prediction = function(x, modeldomain, type, samplesize,variables){

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
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))
    x_clean <- x[complete.cases(x),]

    modeldomain <- modeldomain[,variables]
    modeldomain <- sf::st_drop_geometry(modeldomain)
    modeldomain <- data.frame(scale(modeldomain,center=scaleparam$`scaled:center`,
                                    scale=scaleparam$`scaled:scale`))

    target_dist_feature <- c()
    for (i in 1:nrow(modeldomain)){

      trainDist <-  FNN::knnx.dist(modeldomain[i,],x_clean,k=1)
      target_dist_feature <- c(target_dist_feature,min(trainDist,na.rm=T))
    }
    sampletoprediction <- data.frame(dist = target_dist_feature,
                                     what = "prediction-to-sample",
                                     dist_type = "feature")
  }

  return(sampletoprediction)
}


# sample to test


sample2test <- function(x, testdata, type,variables){

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
    scaleparam <- attributes(scale(x))
    x <- data.frame(scale(x))
    x_clean <- x[complete.cases(x),]
    testdata <- testdata[,variables]
    testdata <- sf::st_drop_geometry(testdata)
    testdata <- data.frame(scale(testdata,center=scaleparam$`scaled:center`,
                                 scale=scaleparam$`scaled:scale`))


    test_dist_feature <- c()
    for (i in 1:nrow(testdata)){

      testDist <- FNN::knnx.dist(testdata[i,],x_clean,k=1)
      test_dist_feature <- c(test_dist_feature,min(testDist,na.rm=T))
    }
    dists_test <- data.frame(dist = test_dist_feature,
                             what = "test-to-sample",
                             dist_type = "feature")


  }
  return(dists_test)
}



# between folds

cvdistance <- function(x, cvfolds, cvtrain, type, variables){

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
    x <- data.frame(scale(x))


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

        trainDist <-  tryCatch(FNN::knnx.dist(testdata_i[k,],traindata_i,k=1),
                               error = function(e)e)
        if(inherits(trainDist, "error")){
          trainDist <- NA
          message("warning: no distance could be calculated for a fold.
                  Possibly because predictor values are NA")
          }

        trainDist[k] <- NA
        d_cv <- c(d_cv,min(trainDist,na.rm=T))
      }
    }


    dists_cv <- data.frame(dist = d_cv,
                           what = factor("CV-distances"),
                           dist_type = "feature")

  }
  return(dists_cv)
}





sampleFromArea <- function(modeldomain, samplesize, type,variables,sampling){

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
    modeldomain <- terra::project(modeldomain, "epsg:4326")
    predictionloc <- sf::st_as_sf(terra::extract(modeldomain,terra::vect(predictionloc),bind=TRUE))
    predictionloc <- na.omit(predictionloc)
  }

  return(predictionloc)

}

# plot results

.plot.nnd <- function(x,type,unit,stat){
  if(unit=="km"){
    x$dist <- x$dist/1000
    xlabs <- "geographic distances (km)"
  }else{
    xlabs <- "geographic distances (m)"
  }

  if( type=="feature"){ xlabs <- "feature space distances"}
  what <- "" #just to avoid check note
  if (type=="feature"){unit ="unitless"}
  if(stat=="density"){
    p <- ggplot2::ggplot(data=x, aes(x=dist, group=what, fill=what)) +
      ggplot2::geom_density(adjust=1.5, alpha=.4, stat=stat) +
      ggplot2::scale_fill_discrete(name = "distance function") +
      ggplot2::xlab(xlabs) +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }else if(stat=="ecdf"){
    p <- ggplot2::ggplot(data=x, aes(x=dist, group=what, col=what)) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step", lwd = 1) +
      ggplot2::scale_color_discrete(name = "distance function") +
      ggplot2::xlab(xlabs) +
      ggplot2::ylab("ECDF") +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }
}
