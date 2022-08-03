#' 'caret' Applications for Spatio-Temporal models
#' @description Supporting functionality to run 'caret' with spatial or spatial-temporal data.
#' 'caret' is a frequently used package for model training and prediction using machine learning.
#' CAST includes functions to improve spatial-temporal modelling tasks using 'caret'.
#' It supports Leave-Location-Out and Leave-Time-Out cross-validation of spatial and spatial-temporal models
#' and allows for spatial variable selection to selects suitable predictor variables
#' in view to their contribution to the spatial model performance.
#' CAST further includes functionality to estimate the (spatial) area of applicability of prediction models
#' by analysing the similarity between new data and training data.
#'
#' @name CAST
#' @docType package
#' @title 'caret' Applications for Spatial-Temporal Models
#' @author Hanna Meyer, Carles Milà, Marvin Ludwig
#' @references
#' \itemize{
#' \item Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
#' \item Meyer, H., Pebesma, E. (2022): Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications. 13.
#' \item Meyer, H., Pebesma, E. (2021): Predicting into unknown space? Estimating the area of applicability of spatial prediction models. Methods in Ecology and Evolution. 12, 1620– 1633.
#' \item Meyer, H., Reudenbach, C., Wöllauer, S., Nauss, T. (2019): Importance of spatial predictor variable selection in machine learning applications - Moving from data reproduction to spatial prediction. Ecological Modelling. 411, 108815.
#' \item Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018): Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software 101: 1-9.
#' }
#'
#' @import caret
#' @importFrom stats sd dist na.omit lm predict quantile na.exclude complete.cases
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @importFrom graphics axis plot segments
#' @keywords package
#'
NULL
