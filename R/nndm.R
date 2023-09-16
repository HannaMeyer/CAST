#' Nearest Neighbour Distance Matching (NNDM) algorithm
#' @description
#' This function implements the NNDM algorithm and returns the necessary
#' indices to perform a NNDM LOO CV for map validation.
#' @author Carles Milà
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param modeldomain sf polygon object defining the prediction area (see Details).
#' @param ppoints sf or sfc point object. Contains the target prediction points.
#' Optional. Alternative to modeldomain (see Details).
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of ppoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See `sf::st_sample`.
#' Only required if modeldomain is used instead of ppoints.
#' @param phi Numeric. Estimate of the landscape autocorrelation range in the
#' same units as the tpoints and ppoints for projected CRS, in meters for geographic CRS.
#' Per default (phi="max"), the size of the prediction area is used. See Details.
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
#' caret's \code{\link{trainControl}} function or used to initiate a custom validation strategy in mlr3.
#'
#' @details NNDM proposes a LOO CV scheme such that the nearest neighbour distance distribution function between the test and training data during the CV process is matched to the nearest neighbour
#' distance distribution function between the prediction and training points. Details of the method can be found in Milà et al. (2022).
#'
#' Specifying \emph{phi} allows limiting distance matching to the area where this is assumed to be relevant due to spatial autocorrelation.
#' Distances are only matched up to \emph{phi}. Beyond that range, all data points are used for training, without exclusions.
#' When \emph{phi} is set to "max", nearest neighbor distance matching is performed for the entire prediction area. Euclidean distances are used for projected
#' and non-defined CRS, great circle distances are used for geographic CRS (units in meters).
#'
#' The \emph{modeldomain} is a sf polygon that defines the prediction area. The function takes a regular point sample (amount defined by \emph{samplesize)} from the spatial extent.
#' As an alternative use \emph{ppoints} instead of \emph{modeldomain}, if you have already defined the prediction locations (e.g. raster pixel centroids).
#' When using either \emph{modeldomain} or \emph{ppoints}, we advise to plot the study area polygon and the training/prediction points as a previous step to ensure they are aligned.
#'
#' @note NNDM is a variation of LOOCV and therefore may take a long time for large training data sets.
#' A k-fold variant will be implemented shortly.
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
#' nndm_pred <- nndm(train_points, ppoints=pred_points)
#' nndm_pred
#' plot(nndm_pred)
#'
#' # ...or run NNDM with a known autocorrelation range of 10
#' # to restrict the matching to distances lower than that.
#' nndm_pred <- nndm(train_points, ppoints=pred_points, phi = 10)
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
#' nndm_pred <- nndm(train_points, ppoints=pred_points)
#' nndm_pred
#' plot(nndm_pred)
#'
#' ########################################################################
#' # Example 3: Real- world example; using a modeldomain instead of previously
#' # sampled prediction locations
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' ### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
#'    by=list(as.character(dat$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#' studyArea[!is.na(studyArea)] <- 1
#' studyArea <- as.polygons(studyArea, values = FALSE, na.all = TRUE) |>
#'     st_as_sf() |>
#'     st_union()
#' pts <- st_transform(pts, crs = st_crs(studyArea))
#' plot(studyArea)
#' plot(st_geometry(pts), add = TRUE, col = "red")
#'
#' nndm_folds <- nndm(pts, modeldomain= studyArea)
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


nndm <- function(tpoints, modeldomain = NULL, ppoints = NULL, samplesize = 1000, sampling = "regular",
                 phi = "max", min_train = 0.5){

  # create sample points from modeldomain
  if(is.null(ppoints)&!is.null(modeldomain)){
    if(!identical(sf::st_crs(tpoints), sf::st_crs(modeldomain))){
      stop("tpoints and modeldomain must have the same CRS")
    }
    message(paste0(samplesize, " prediction points are sampled from the modeldomain"))
    ppoints <- sf::st_sample(x = modeldomain, size = samplesize, type = sampling)
    sf::st_crs(ppoints) <- sf::st_crs(modeldomain)
  }else if(!is.null(ppoints)){
    if(!identical(sf::st_crs(tpoints), sf::st_crs(ppoints))){
      stop("tpoints and ppoints must have the same CRS")
    }
  }

  # If tpoints is sfc, coerce to sf.
  if(any(class(tpoints) %in% "sfc")){
    tpoints <- sf::st_sf(geom=tpoints)
  }

  # If ppoints is sfc, coerce to sf.
  if(any(class(ppoints) %in% "sfc")){
    ppoints <- sf::st_sf(geom=ppoints)
  }

  # if phi==max calculate the range of the size area
  if(phi=="max"){
    xmin <- min(sf::st_coordinates(ppoints)[,1])
    xmax <- max(sf::st_coordinates(ppoints)[,1])
    ymin <- min(sf::st_coordinates(ppoints)[,2])
    ymax <-  max(sf::st_coordinates(ppoints)[,2])
    p <- sf::st_sfc(sf::st_point(c(xmin,ymin)), sf::st_point(c(xmax,ymax)))
    sf::st_crs(p) <- sf::st_crs(ppoints)
    phi <- as.numeric(max(sf::st_distance(p)))
  }

  # Input data checks
  nndm_checks(tpoints, ppoints, phi, min_train)

  # Compute nearest neighbour distances between training and prediction points
  Gij <- sf::st_distance(ppoints, tpoints)
  units(Gij) <- NULL
  Gij <- apply(Gij, 1, min)

  # Compute distance matrix of training points
  tdist <- sf::st_distance(tpoints)
  units(tdist) <- NULL
  diag(tdist) <- NA
  Gj <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
  Gjstar <- Gj

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
nndm_checks <- function(tpoints, ppoints, phi, min_train){

  # Check for valid range of phi
  if(phi < 0 | !is.numeric(phi)){
    stop("phi must be positive.")
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

  # Check class and geometry type of ppoints
  if(!any(c("sfc", "sf") %in% class(ppoints))){
    stop("ppoints must be a sf/sfc object.")
  }else if(!any(class(sf::st_geometry(ppoints)) %in% c("sfc_POINT"))){
    stop("ppoints must be a sf/sfc point object.")
  }

}
