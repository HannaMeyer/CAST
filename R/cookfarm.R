#' Cookfarm soil logger data
#'
#' spatio-temporal data of soil properties and associated predictors for the Cookfarm in Washington, USA
#' @format
#' A sf data.frame with 128545 rows and 17 columns:
#' \describe{
#'   \item{SOURCEID}{ID of the logger}
#'   \item{VW}{Response Variable - Soil Moisture}
#'   \item{altitude}{Measurement depth of VW}
#'   \item{Date, cdata}{Measurement Date, Cumulative Date}
#'   \item{Easting, Northing}{Location Coordinates (EPSG:26911)}
#'   \item{DEM, TWI, NDRE.M, NDRE.Sd, Precip_wrcc, MaxT_wrcc, MinT_wrcc, Precip_cum}{Predictor Variables}
#' }
#'
#' @references \itemize{
#' \item{Gash et al. 2015 - Spatio-temporal interpolation of soil water, temperature, and electrical conductivity in 3D + T: The Cook Agronomy Farm data set \doi{https://doi.org/10.1016/j.spasta.2015.04.001}}
#' \item{Meyer et al. 2018 - Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation \doi{https://doi.org/10.1016/j.envsoft.2017.12.001}}
#' }
#' @usage data(cookfarm)
#'
"cookfarm"

