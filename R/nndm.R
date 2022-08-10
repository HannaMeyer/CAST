#' Nearest Neighbour Distance Matching (NNDM) algorithm
#' @description
#' This function implements the \emph{NNDM} algorithm and returns the necessary
#' indices to perform a NNDM LOO CV for map validation.
#' @author Carles Milà
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param ppoints sf or sfc point object. Contains the target prediction points.
#' @param phi Numeric. Estimate of the landscape autocorrelation range in the
#' same units as the tpoints and ppoints for projected CRS, in meters for geographic CRS.
#' Per default (phi="max"), the size of the prediction area is used. See Details
#' @param min_train Numeric between 0 and 1. Minimum proportion of training
#' data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions).
#'
#' @return An object of class \emph{nndm} consisting of a list of six elements:
#' indx_train, indx_test, and indx_exclude (indices of the observations to use as
#' training/test/excluded data in each NNDM LOO CV iteration), Gij (distances for
#' multitype G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during NNDM LOO CV), phi (landscape autocorrelation range).
#' indx_train and indx_test can directly be used as "index" and "indexOut" in
#' caret's \code{\link{trainControl}} function.
#' @details Details of the method can be found in Milà et al. (2022).
#' Euclidean distances are used for projected
#' and non-defined CRS, great circle distances are used for geographic CRS (units in meters).
#' Specifying phi allows limiting distance matching to the area where this is assumed to be relevant due to spatial autocorrelation.
#' Distances are only matched up to phi. Beyond that range, all data points are used for training, without exclusions.
#' When phi is set to "max", nearest neighbor distance matching is performed for the entire prediction area.
#' @note NNDM is a variation of LOOCV and therefore may take a long time for large training data sets.
#' You may need to consider alternatives following the ideas of Milà et al. (2022) for large data sets.
#' @seealso \code{\link{plot_geodist}}
#' @references
#' \itemize{
#' \item Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
#' \item Meyer, H., Pebesma, E. (2022): Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications. 13.
#' }
#' @export
#' @examples
#' library(sf)
#'
#' # Simulate 100 random training and test points in a 100x100 square
#' set.seed(123)
#' poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' sample_poly <- sf::st_polygon(poly)
#' train_points <- sf::st_sample(sample_poly, 100, type = "random")
#' pred_points <- sf::st_sample(sample_poly, 100, type = "random")
#'
#' # Run NNDM.
#' nndm_pred <- nndm(train_points, pred_points)
#' nndm_pred
#' plot(nndm_pred)
#'
#' # ...or run NNDM with a known autocorrelation range.
#' # Here, the autocorrelation range (phi) is known to be 10.
#' nndm_pred <- nndm(train_points, pred_points, 10)
#' nndm_pred
#' plot(nndm_pred)

nndm <- function(tpoints, ppoints, phi="max", min_train=0){

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

  # min_train must be a single positive integer
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

  # Check same CRS of tpoints and ppoints
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
    stop("tpoints and ppoints must have the same CRS.")
  }
}
