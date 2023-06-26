#' sPlotOpen Data of Species Richness
#'
#' sPlotOpen Species Richness for South America with associated predictors from worldclim
#' @format
#' A sf points / data.frame with 703 rows and 17 columns:
#' \describe{
#'   \item{PlotObeservationID, GIVD_ID, Country, Biome}{sPlotOpen Metadata}
#'   \item{Species_richness}{Response Variable - Plant species richness from sPlotOpen}
#'   \item{bio_x, elev}{Predictor Variables - Worldclim and SRTM}
#'   \item{geometry}{Lat/Lon}
#' }
#' @source sPlotOpen: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13346
#'     predictors acquired via R package geodata
#' @usage data(splotdata)
#'
"splotdata"
