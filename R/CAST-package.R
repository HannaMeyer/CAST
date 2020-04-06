#' 'caret' Applications for Spatio-Temporal models
#'
#' @description Supporting functionality to run 'caret' with spatial or spatial-temporal data.
#' 'caret' is a frequently used package for model training and prediction using machine learning.
#' CAST includes functions to improve spatial-temporal modelling tasks using 'caret'.
#' It supports Leave-Location-Out and Leave-Time-Out cross-validation of spatial and spatial-temporal models
#' and allows for spatial variable selection to selects suitable predictor variables
#' in view to their contribution to the spatial model performance.
#' CAST further includes functionality to estimate the (spatial) area of applicability of prediction models
#' by analysing the similarity between new data and training data.

#' @name CAST
#' @docType package
#' @title 'caret' Applications for Spatial-Temporal Models
#' @author Hanna Meyer, Christoph Reudenbach, Marvin Ludwig, Thomas Nauss
#' \cr
#' \emph{Maintainer:} Hanna Meyer \email{hanna.meyer@@uni-muenster.de}
#'
#' @import caret
#' @importFrom stats sd dist na.omit lm predict quantile na.exclude
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom graphics axis plot segments
#' @keywords package
#'
NULL
