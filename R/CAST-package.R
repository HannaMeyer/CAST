#' 'caret' Applications for Spatio-Temporal models
#'
#' @description Supporting functionality to run 'caret' with spatial or spatial-temporal data. 'caret' is a frequently used package for model training and prediction using machine learning. This package includes functions to improve spatial-temporal modelling tasks using 'caret'. It prepares data for Leave-Location-Out and Leave-Time-Out cross-validation which are target-oriented validation strategies for spatial-temporal models. To decrease overfitting and improve model performances, the package implements a forward feature selection that selects suitable predictor variables in view to their contribution to the target-oriented performance.

#' @name CAST
#' @docType package
#' @title 'caret' Applications for Spatial-Temporal Models
#' @author Hanna Meyer, Christoph Reudenbach, Marvin Ludwig, Thomas Nauss
#' \cr
#' \emph{Maintainer:} Hanna Meyer \email{hanna.meyer@@uni-muenster.de}
#'
#' @import caret

#' @importFrom stats na.exclude
#' @importFrom stats sd dist na.omit
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom graphics axis plot segments
#' @keywords package
#'
NULL
