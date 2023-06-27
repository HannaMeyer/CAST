#' sPlotOpen Data of Species Richness
#'
#' sPlotOpen Species Richness for South America with associated predictors
#' @format
#' A sf points / data.frame with 703 rows and 17 columns:
#' \describe{
#'   \item{PlotObeservationID, GIVD_ID, Country, Biome}{sPlotOpen Metadata}
#'   \item{Species_richness}{Response Variable - Plant species richness from sPlotOpen}
#'   \item{bio_x, elev}{Predictor Variables - Worldclim and SRTM elevation}
#'   \item{geometry}{Lat/Lon}
#' }
#' @source \itemize{
#' \item{Plot with Species_richness from sPlotOpen: \link{https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13346}}
#' \item{predictors acquired via R package geodata \link{https://github.com/rspatial/geodata}}
#' }
#'
#' @references \itemize{
#' \item{Sabatini, F. M. et al. sPlotOpen – An environmentally balanced, open‐access, global dataset of vegetation plots. (2021). \doi{10.1111/geb.13346}}
#' \item{Lopez-Gonzalez, G. et al. ForestPlots.net: a web application and research tool to manage and analyse tropical forest plot data: ForestPlots.net.
#'  Journal of Vegetation Science (2011).}
#' \item{Pauchard, A. et al. Alien Plants Homogenise Protected Areas: Evidence from the Landscape and Regional Scales in South Central Chile. in Plant Invasions in Protected Areas (2013).}
#' \item{Peyre, G. et al. VegPáramo, a flora and vegetation database for the Andean páramo. phytocoenologia (2015).}
#' \item{Vibrans, A. C. et al. Insights from a large-scale inventory in the southern Brazilian Atlantic Forest. Scientia Agricola (2020).}
#' }
#' @usage data(splotdata)
#'
"splotdata"
