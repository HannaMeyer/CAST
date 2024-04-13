#' 'caret' Applications for Spatio-Temporal models
#' @description Supporting functionality to run 'caret' with spatial or spatial-temporal data.
#' 'caret' is a frequently used package for model training and prediction using machine learning.
#' CAST includes functions to improve spatial-temporal modelling tasks using 'caret'.
#' It includes the newly suggested 'Nearest neighbor distance matching' cross-validation to estimate the performance
#' of spatial prediction models and allows for spatial variable selection to selects suitable predictor variables
#' in view to their contribution to the spatial model performance.
#' CAST further includes functionality to estimate the (spatial) area of applicability of prediction models
#' by analysing the similarity between new data and training data.
#' Methods are described in Meyer et al. (2018); Meyer et al. (2019); Meyer and Pebesma (2021); Milà et al. (2022); Meyer and Pebesma (2022); Linnenbrink et al. (2023).
#' The package is described in detail in Meyer et al. (2024).
#' @name CAST
#' @title 'caret' Applications for Spatial-Temporal Models
#' @author Hanna Meyer, Carles Milà, Marvin Ludwig, Jan Linnenbrink, Fabian Schumacher
#' @references
#' \itemize{
#' \item Meyer, H., Ludwig, L., Milà, C., Linnenbrink, J., Schumacher, F. (2024): The CAST package for training and assessment of spatial prediction models in R. arXiv, https://doi.org/10.48550/arXiv.2404.06978.
#' \item Linnenbrink, J., Milà, C., Ludwig, M., and Meyer, H.: kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1308, 2023.
#' \item Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
#' \item Meyer, H., Pebesma, E. (2022): Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications. 13.
#' \item Meyer, H., Pebesma, E. (2021): Predicting into unknown space? Estimating the area of applicability of spatial prediction models. Methods in Ecology and Evolution. 12, 1620– 1633.
#' \item Meyer, H., Reudenbach, C., Wöllauer, S., Nauss, T. (2019): Importance of spatial predictor variable selection in machine learning applications - Moving from data reproduction to spatial prediction. Ecological Modelling. 411, 108815.
#' \item Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018): Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software 101: 1-9.
#' }
#'
#' @import caret
#' @importFrom stats sd dist na.omit lm predict quantile na.exclude complete.cases median
#' @importFrom utils combn txtProgressBar setTxtProgressBar
#' @importFrom grDevices rainbow
#' @importFrom graphics axis plot segments
#' @keywords package
#' @aliases CAST-package
#'
"_PACKAGE"
