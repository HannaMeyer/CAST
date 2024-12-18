#' Nearest Neighbour Distance Matching (NNDM) algorithm
#' @description
#' This function implements the NNDM algorithm and returns the necessary indices to perform a NNDM LOO CV for map validation.
#' @author Carles Milà
#' @param tpoints sf or sfc point object, or data.frame if space = "feature". Contains the training points samples.
#' @param modeldomain sf polygon object or SpatRaster defining the prediction area. Optional; alternative to predpoints (see Details).
#' @param predpoints sf or sfc point object, or data.frame if space = "feature". Contains the target prediction points. Optional; alternative to modeldomain (see Details).
#' @param space character. Either "geographical" or "feature". Feature space is still experimental, so use with caution.
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of predpoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See `sf::st_sample`.
#' Only required if modeldomain is used instead of predpoints.
#' @param phi Numeric. Estimate of the landscape autocorrelation range in the
#' same units as the tpoints and predpoints for projected CRS, in meters for geographic CRS.
#' Per default (phi="max"), the maximum distance found in the training and prediction points is used. See Details.
#' @param min_train Numeric between 0 and 1. Minimum proportion of training
#' data that must be used in each CV fold. Defaults to 0.5 (i.e. half of the training points).
#'
#' @return An object of class \emph{nndm} consisting of a list of six elements:
#' indx_train, indx_test, and indx_exclude (indices of the observations to use as
#' training/test/excluded data in each NNDM LOO CV iteration), Gij (distances for
#' G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during NNDM LOO CV), phi (landscape autocorrelation range).
#' indx_train and indx_test can directly be used as "index" and "indexOut" in
#' caret's \code{\link[caret]{trainControl}} function or used to initiate a custom validation strategy in mlr3.
#'
#' @details NNDM proposes a LOO CV scheme such that the nearest neighbour distance distribution function between the test and training data during the CV process is matched to the nearest neighbour
#' distance distribution function between the prediction and training points. Details of the method can be found in Milà et al. (2022).
#'
#' Specifying \emph{phi} allows limiting distance matching to the area where this is assumed to be relevant due to spatial autocorrelation.
#' Distances are only matched up to \emph{phi}. Beyond that range, all data points are used for training, without exclusions.
#' When \emph{phi} is set to "max", nearest neighbor distance matching is performed for the entire prediction area. Euclidean distances are used for projected
#' and non-defined CRS, great circle distances are used for geographic CRS (units in meters).
#'
#' The \emph{modeldomain} is either a sf polygon that defines the prediction area, or alternatively a SpatRaster out of which a polygon,
#' transformed into the CRS of the training points, is defined as the outline of all non-NA cells.
#' Then, the function takes a regular point sample (amount defined by \emph{samplesize)} from the spatial extent.
#' As an alternative use \emph{predpoints} instead of \emph{modeldomain}, if you have already defined the prediction locations (e.g. raster pixel centroids).
#' When using either \emph{modeldomain} or \emph{predpoints}, we advise to plot the study area polygon and the training/prediction points as a previous step to ensure they are aligned.
#'
#' @note NNDM is a variation of LOOCV and therefore may take a long time for large training data sets. See \code{\link{knndm}} for a more efficient k-fold variant of the method.
#' @seealso \code{\link{geodist}}, \code{\link{knndm}}
#' @references
#' \itemize{
#' \item Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
#' \item Meyer, H., Pebesma, E. (2022): Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications. 13.
#' }
#' @export
#' @examples
#' ########################################################################
#' # Example 1: Simulated data - Randomly-distributed training points
#' ########################################################################
#'
#' library(sf)
#'
#' # Simulate 100 random training points in a 100x100 square
#' set.seed(123)
#' poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' sample_poly <- sf::st_polygon(poly)
#' train_points <- sf::st_sample(sample_poly, 100, type = "random")
#' pred_points <- sf::st_sample(sample_poly, 100, type = "regular")
#' plot(sample_poly)
#' plot(pred_points, add = TRUE, col = "blue")
#' plot(train_points, add = TRUE, col = "red")
#'
#' # Run NNDM for the whole domain, here the prediction points are known
#' nndm_pred <- nndm(train_points, predpoints=pred_points)
#' nndm_pred
#' plot(nndm_pred)
#' plot(nndm_pred, type = "simple") # For more accessible legend labels
#'
#' # ...or run NNDM with a known autocorrelation range of 10
#' # to restrict the matching to distances lower than that.
#' nndm_pred <- nndm(train_points, predpoints=pred_points, phi = 10)
#' nndm_pred
#' plot(nndm_pred)
#'
#' ########################################################################
#' # Example 2: Simulated data - Clustered training points
#' ########################################################################
#'
#' library(sf)
#'
#' # Simulate 100 clustered training points in a 100x100 square
#' set.seed(123)
#' poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' sample_poly <- sf::st_polygon(poly)
#' train_points <- clustered_sample(sample_poly, 100, 10, 5)
#' pred_points <- sf::st_sample(sample_poly, 100, type = "regular")
#' plot(sample_poly)
#' plot(pred_points, add = TRUE, col = "blue")
#' plot(train_points, add = TRUE, col = "red")
#'
#' # Run NNDM for the whole domain
#' nndm_pred <- nndm(train_points, predpoints=pred_points)
#' nndm_pred
#' plot(nndm_pred)
#' plot(nndm_pred, type = "simple") # For more accessible legend labels
#'
#' ########################################################################
#' # Example 3: Real- world example; using a SpatRast modeldomain instead
#' # of previously sampled prediction locations
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' ### prepare sample data:
#' data(cookfarm)
#' dat <- aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
#'    by=list(as.character(cookfarm$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
#' pts <- st_transform(pts, crs = st_crs(studyArea))
#' terra::plot(studyArea[["DEM"]])
#' terra::plot(vect(pts), add = T)
#'
#' nndm_folds <- nndm(pts, modeldomain = studyArea)
#' plot(nndm_folds)
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=nndm_folds$indx_train,
#'    indexOut=nndm_folds$indx_test,
#'    savePredictions='final')
#' model_nndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
#'    dat$VW,
#'    method="rf",
#'    trControl = ctrl)
#' global_validation(model_nndm)
#'}
#'
#' ########################################################################
#' # Example 4: Real- world example; nndm in feature space
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#'
#' # Prepare the splot dataset for Chile
#' data(splotdata)
#' splotdata <- splotdata[splotdata$Country == "Chile",]
#'
#' # Select a series of bioclimatic predictors
#' predictors <- c("bio_1", "bio_4", "bio_5", "bio_6",
#'                "bio_8", "bio_9", "bio_12", "bio_13",
#'                "bio_14", "bio_15", "elev")
#'
#' predictors_sp <- terra::rast(system.file("extdata", "predictors_chile.tif", package="CAST"))
#'
#' # Data visualization
#' terra::plot(predictors_sp[["bio_1"]])
#' terra::plot(vect(splotdata), add = T)
#'
#' # Run and visualise the nndm results
#' nndm_folds <- nndm(splotdata[,predictors], modeldomain = predictors_sp, space = "feature")
#' plot(nndm_folds)
#'
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=nndm_folds$indx_train,
#'    indexOut=nndm_folds$indx_test,
#'    savePredictions='final')
#' model_nndm <- train(st_drop_geometry(splotdata[,predictors]),
#'    splotdata$Species_richness,
#'    method="rf",
#'    trControl = ctrl)
#' global_validation(model_nndm)
#'
#' }
nndm <- function(tpoints, modeldomain = NULL, predpoints = NULL,
                 space="geographical",
                 samplesize = 1000, sampling = "regular",
                 phi = "max", min_train = 0.5){


  # 1. Preprocessing actions ----
  if(is.null(predpoints)&!is.null(modeldomain)){

    # Check modeldomain is indeed a sf/SpatRaster
    if(!any(c("sfc", "sf", "SpatRaster") %in% class(modeldomain))){
      stop("modeldomain must be a sf/sfc object or a 'SpatRaster' object.")
    }

    # If modeldomain is a SpatRaster, transform into polygon
    if(any(class(modeldomain) == "SpatRaster")){

      # save predictor stack for extraction if space = "feature"
      if(space == "feature") {
        predictor_stack <- modeldomain
      }
      modeldomain[!is.na(modeldomain)] <- 1
      modeldomain <- terra::as.polygons(modeldomain, values = FALSE, na.all = TRUE) |>
        sf::st_as_sf() |>
        sf::st_union()
      if(any(c("sfc", "sf") %in% class(tpoints))) {
        modeldomain <- sf::st_transform(modeldomain, crs = sf::st_crs(tpoints))
      }
    }

    # Check modeldomain is indeed a polygon sf
    if(!any(class(sf::st_geometry(modeldomain)) %in% c("sfc_POLYGON", "sfc_MULTIPOLYGON"))){
      stop("modeldomain must be a sf/sfc polygon object.")
    }

    # Check whether modeldomain has the same crs as tpoints
    if(!identical(sf::st_crs(tpoints), sf::st_crs(modeldomain)) & space == "geographical"){
      stop("tpoints and modeldomain must have the same CRS")
    }

    # We sample
    message(paste0(samplesize, " prediction points are sampled from the modeldomain"))
    predpoints <- sf::st_sample(x = modeldomain, size = samplesize, type = sampling)
    sf::st_crs(predpoints) <- sf::st_crs(modeldomain)

    if(space == "feature") {
      message("predictor values are extracted for prediction points")
      predpoints <- terra::extract(predictor_stack, terra::vect(predpoints), ID=FALSE)
    }

  }else if(!is.null(predpoints) & space == "geographical"){
    if(!identical(sf::st_crs(tpoints), sf::st_crs(predpoints))){
      stop("tpoints and predpoints must have the same CRS")
    }
  }


  # 2. Data formats, data checks, scaling and categorical variables ----
  if(isTRUE(space == "geographical")) {

    # If tpoints is sfc, coerce to sf.
    if(any(class(tpoints) %in% "sfc")){
      tpoints <- sf::st_sf(geom=tpoints)
    }
    # If predpoints is sfc, coerce to sf.
    if(any(class(predpoints) %in% "sfc")){
      predpoints <- sf::st_sf(geom=predpoints)
    }

    # Input data checks
    nndm_checks_geo(tpoints, predpoints, phi, min_train)

  }else if(isTRUE(space == "feature")){

    # drop geometry if tpoints / predpoints are of class sf
    if(any(class(tpoints) %in% c("sf","sfc"))) {
      tpoints <- sf::st_set_geometry(tpoints, NULL)
    }
    if(any(class(predpoints) %in% c("sf","sfc"))) {
      predpoints <- sf::st_set_geometry(predpoints, NULL)
    }

    # get names of categorical variables
    catVars <- names(tpoints)[which(sapply(tpoints, class)%in%c("factor","character"))]
    if(length(catVars)==0) {
      catVars <- NULL
    }
    if(!is.null(catVars)) {
      message(paste0("variable(s) '", catVars, "' is (are) treated as categorical variables"))
    }

    # omit NAs
    if(any(is.na(predpoints))) {
      message("some prediction points contain NAs, which will be removed")
      predpoints <- stats::na.omit(predpoints)
    }
    if(any(is.na(tpoints))) {
      message("some training points contain NAs, which will be removed")
      tpoints <- stats::na.omit(tpoints)
    }

    # Input data checks
    nndm_checks_feature(tpoints, predpoints, phi, min_train, catVars)

    # Scaling and dealing with categorical factors
    if(is.null(catVars)) {

      scale_attr <- attributes(scale(tpoints))
      tpoints <- scale(tpoints) |> as.data.frame()
      predpoints <- scale(predpoints,center=scale_attr$`scaled:center`,
                          scale=scale_attr$`scaled:scale`) |>
        as.data.frame()

    } else {
      tpoints_cat <- tpoints[,catVars,drop=FALSE]
      predpoints_cat <- predpoints[,catVars,drop=FALSE]

      tpoints_num <- tpoints[,-which(names(tpoints)%in%catVars),drop=FALSE]
      predpoints_num <- predpoints[,-which(names(predpoints)%in%catVars),drop=FALSE]

      scale_attr <- attributes(scale(tpoints_num))
      tpoints <- scale(tpoints_num) |> as.data.frame()
      predpoints <- scale(predpoints_num,center=scale_attr$`scaled:center`,
                          scale=scale_attr$`scaled:scale`) |>
        as.data.frame()
      tpoints <- as.data.frame(cbind(tpoints, lapply(tpoints_cat, as.factor)))
      predpoints <- as.data.frame(cbind(predpoints, lapply(predpoints_cat, as.factor)))

      # 0/1 encode categorical variables (as in R/trainDI.R)
      for (catvar in catVars){
        # mask all unknown levels in newdata as NA
        tpoints[,catvar]<-droplevels(tpoints[,catvar])
        predpoints[,catvar]<-droplevels(predpoints[,catvar])

        # then create dummy variables for the remaining levels in train:
        dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = tpoints),
                             tpoints)
        dvi_predpoints <- predict(caret::dummyVars(paste0("~",catvar), data = predpoints),
                                  predpoints)
        tpoints <- data.frame(tpoints,dvi_train)
        predpoints <- data.frame(predpoints,dvi_predpoints)

      }
      tpoints <- tpoints[,-which(names(tpoints)%in%catVars)]
      predpoints <- predpoints[,-which(names(predpoints)%in%catVars)]

    }
  }

  # 3. Distance and phi computation ----
  if(isTRUE(space=="geographical")){

    # Compute nearest neighbour distances between training and prediction points
    Gij <- sf::st_distance(predpoints, tpoints)
    units(Gij) <- NULL
    Gij <- apply(Gij, 1, min)

    # Compute distance matrix of training points
    tdist <- sf::st_distance(tpoints)
    units(tdist) <- NULL
    diag(tdist) <- NA
    Gj <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
    Gjstar <- Gj

    # if phi==max calculate the maximum relevant distance
    if(phi=="max"){
      phi <- max(c(Gij, c(tdist)), na.rm=TRUE) + 1e-9
    }


  }else if(isTRUE(space=="feature")){

    if(is.null(catVars)) {

      # Euclidean distances if no categorical variables are present
      Gij <- c(FNN::knnx.dist(query = predpoints, data = tpoints, k = 1))
      tdist <- as.matrix(stats::dist(tpoints, upper = TRUE))
      diag(tdist) <- NA
      Gj <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
      Gjstar <- Gj

    } else {

      # Gower distances if categorical variables are present
      Gj <- sapply(1:nrow(tpoints), function(i) gower::gower_topn(tpoints[i,], tpoints[-i,], n=1)$distance[[1]])
      tdist <- matrix(NA, nrow=nrow(tpoints), ncol=nrow(tpoints))
      for(r in 1:nrow(tdist)){
        tdist[r,] <- gower::gower_dist(tpoints[r,], tpoints)
      }
      diag(tdist) <- NA
      Gj <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
      Gjstar <- Gj

    }

    # if phi==max calculate the maximum relevant distance
    if(phi=="max"){
      phi <- max(c(Gij, c(tdist)), na.rm=TRUE) + 1e-9
    }
  }

  # Start algorithm
  rmin <- min(Gjstar)
  jmin <- which.min(Gjstar)[1]
  kmin <- which(tdist[jmin,]==rmin)

  while(rmin <= phi){

    # Check if removing the point improves the match. If yes, update
    if((sum(Gjstar<=rmin)-1)/length(Gjstar) >= (sum(Gij<=rmin)/length(Gij)) &
       sum(!is.na(tdist[jmin, ]))/ncol(tdist) > min_train){
      tdist[jmin, kmin] <- NA
      Gjstar <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
      rmin <- min(Gjstar[Gjstar>=rmin]) # Distances are the same for the same pair
      jmin <- which(Gjstar==rmin)[1]
      kmin <- which(tdist[jmin,]==rmin)

    }else if(sum(Gjstar>rmin)==0){
      break
    }else{ # Otherwise move on to the next distance
      rmin <- min(Gjstar[Gjstar>rmin])
      jmin <- which(Gjstar==rmin)[1]
      kmin <- which(tdist[jmin,]==rmin)
    }
  }

  # Derive indicators
  indx_train <- list()
  indx_test <- list()
  indx_exclude <- list()
  for(i in 1:nrow(tdist)){
    indx_train[[i]] <- which(!is.na(tdist[i,]))
    indx_test[[i]] <- i
    indx_exclude[[i]] <- setdiff(which(is.na(tdist[i,])), i)
  }

  # Return list of indices
  res <- list(indx_train=indx_train, indx_test=indx_test,
              indx_exclude=indx_exclude, Gij=Gij, Gj=Gj, Gjstar=Gjstar, phi=phi)
  class(res) <- c("nndm", "list")
  res

}


# Input data checks for NNDM
nndm_checks_geo <- function(tpoints, predpoints, phi, min_train){

  # Check for valid range of phi
  if(phi < 0 | (!is.numeric(phi) & phi!= "max")){
    stop("phi must be positive or set to 'max'.")
  }

  # min_train must be a single positive numeric
  if(length(min_train)!=1 | min_train<0 | min_train>1 | !is.numeric(min_train)){
    stop("min_train must be a numeric between 0 and 1.")
  }

  # Check class and geometry type of tpoints
  if(!any(c("sfc", "sf") %in% class(tpoints))){
    stop("tpoints must be a sf/sfc object.")
  }else if(!any(class(sf::st_geometry(tpoints)) %in% c("sfc_POINT"))){
    stop("tpoints must be a sf/sfc point object.")
  }

  # Check class and geometry type of predpoints
  if(!any(c("sfc", "sf") %in% class(predpoints))){
    stop("predpoints must be a sf/sfc object.")
  }else if(!any(class(sf::st_geometry(predpoints)) %in% c("sfc_POINT"))){
    stop("predpoints must be a sf/sfc point object.")
  }

}

nndm_checks_feature <- function(tpoints, predpoints, phi, min_train, catVars){

  # Check for valid range of phi
  if(phi < 0 | (!is.numeric(phi) & phi!= "max")){
    stop("phi must be positive or set to 'max'.")
  }

  # min_train must be a single positive numeric
  if(length(min_train)!=1 | min_train<0 | min_train>1 | !is.numeric(min_train)){
    stop("min_train must be a numeric between 0 and 1.")
  }

  if(length(setdiff(names(tpoints), names(predpoints)))>0) {
    stop("tpoints and predpoints need to contain the predictor data and have the same colnames.")
  }

  for (catvar in catVars) {
    if (any(!unique(tpoints[,catvar]) %in% unique(predpoints[,catvar]))) {
      stop(paste0("Some values of factor", catvar, "are only present in training / prediction points.
                  All factor values in the prediction points must be present in the training points."))
    }
  }
}
