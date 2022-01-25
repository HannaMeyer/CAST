#' Plot euclidean nearest neighbor distances in geographic space
#'
#' @description Density plot of geographic nearest neighbor distances between training data as well as between training data and prediction locations.
#' Optional, the nearest neighbor distances between training data and test data or between training data and CV iterations is shown.
#' The plot can be used to check the suitability of a chosen CV method to be representative to estimate map accuracy.
#'
#' @param x object of class sf or TrainDI
#' @param modeldomain sf or raster object defining the prediction area
#' @param cvfolds optional. List of row indices of x that are held back in each CV iteration. See e.g. ?createFolds or ?createSpaceTimeFolds
#' @param testdata optional. object of class sf: Data used for independent validation
#' @param samplesize numeric. How many random prediction samples should be used?
#' @param scale logical. Present distances on log scale?
#' @return A plot and a data.frame containing the spatial distances
#'
#' @import ggplot2
#'
#' @author Hanna Meyer, Edzer Pebesma
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#' ########### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' pts_train <- pts[1:29,]
#' pts_test <- pts[30:42,]
#'
#' studyArea <- raster(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#'
#'
#' ########### Distance between training data and new data:
#' dist <- plot_geodist(pts_train,studyArea)
#'
#' ########### Distance between training data, new data and test data:
#' mapview(pts_train,col.regions="blue")+mapview(pts_test,col.regions="red")
#' dist <- plot_geodist(pts_train,studyArea,testdata=pts_test)
#'
#' ########### Distance between training data, new data and CV folds:
#' folds <- createFolds(1:nrow(pts_train),k=3,returnTrain=FALSE)
#' plot_geodist(x=pts,modeldomain=studyArea,cvfolds=folds)
#' }
#'
#'
#'
#'
#'
#' @export


plot_geodist <- function(x,modeldomain,samplesize=2000,
                         cvfolds=NULL,testdata=NULL,scale=TRUE){
  x <- sf::st_transform(x,4326)

  ##### Distance within training data:
  sf::sf_use_s2(TRUE)
  d = sf::st_distance(x)
  diag(d) = Inf
  min_d = apply(d, 1, min)

  ##### Distance to prediction locations:
  # regularly spread points (prediction locations):
  # see https://edzer.github.io/OGH21/

  if (inherits(modeldomain, "Raster")){
    if(samplesize>raster::ncell(modeldomain)){
      samplesize <- raster::ncell(modeldomain)
      message(paste0("samplesize for new data shouldn't be larger than number of pixels.
              Samplesize was reduced to ",raster::ncell(modeldomain)))
    }
    modeldomain <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(modeldomain)))

  }
  sf::sf_use_s2(FALSE)
  sf::st_as_sf(modeldomain) |>
    sf::st_transform(4326) -> bb
  methods::as(bb, "Spatial") |>
    sp::spsample(n =samplesize, type = "regular")  |>
    sf::st_as_sfc() |>
    sf::st_set_crs(4326) -> predictionloc

  sf::sf_use_s2(TRUE)
  d0 <- sf::st_distance(predictionloc, x)
  min_d0 <- apply(d0, 1, min)

  ##### Compile data:
  what <- c()
  dists <- data.frame(dist = c(min_d, min_d0), what = c(rep("sample-to-sample", length(min_d)), rep("sample-to-predition", length(min_d0))))
  dists$what <- factor(dists$what, levels = c("sample-to-sample", "sample-to-predition"))


  ##### Distance to test data:
  if(!is.null(testdata)){
    testdata <- sf::st_transform(testdata,4326)
    d_test <- sf::st_distance(testdata, x)
    min_d_test <- apply(d_test, 1, min)

    dists_test <- data.frame(dist = min_d_test, what = rep("sample-to-test", length(min_d_test)))
    dists_test$what <- factor(dists_test$what, levels = c("sample-to-test"))
    dists = rbind(dists, dists_test)

  }
  ##### Distance to CV data:
  if(!is.null(cvfolds)){
    d_cv <- c()
    for (i in 1:length(cvfolds)){
      d_cv_tmp <- sf::st_distance(x[cvfolds[[i]],], x[-cvfolds[[i]],])
      d_cv <-c(d_cv,apply(d_cv_tmp, 1, min))
    }

    dists_cv <- data.frame(dist = d_cv, what = rep("CV-distances", length(d_cv)))
    dists_cv$what <- factor(dists_cv$what, levels = c("CV-distances"))
    dists <- rbind(dists, dists_cv)

  }

  ##### Compile and plot data:

  p <- ggplot2::ggplot(data=dists, aes(x=dist, group=what, fill=what)) +
    ggplot2::geom_density(adjust=1.5, alpha=.4) +
    ggplot2::scale_fill_discrete(name = "distance function") + xlab("distance (m)") +
    ggplot2::theme(legend.position="bottom")

  if(scale){
    p <- p+scale_x_log10(labels=round)
  }

  print(p)
  return(dists)
}

