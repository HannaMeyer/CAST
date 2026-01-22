# sPlotOpen Data of Species Richness

sPlotOpen Species Richness for South America with associated predictors

## Usage

``` r
data(splotdata)
```

## Format

A sf points / data.frame with 703 rows and 17 columns:

- PlotObeservationID, GIVD_ID, Country, Biome:

  sPlotOpen Metadata

- Species_richness:

  Response Variable - Plant species richness from sPlotOpen

- bio_x, elev:

  Predictor Variables - Worldclim and SRTM elevation

- geometry:

  Lat/Lon

## Source

- Plot with Species_richness from
  [sPlotOpen](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13346)

- predictors acquired via R package
  [geodata](https://github.com/rspatial/geodata)

## References

- Sabatini, F. M. et al. sPlotOpen – An environmentally balanced,
  open‐access, global dataset of vegetation plots. (2021).
  [doi:10.1111/geb.13346](https://doi.org/10.1111/geb.13346)

- Lopez-Gonzalez, G. et al. ForestPlots.net: a web application and
  research tool to manage and analyse tropical forest plot data:
  ForestPlots.net. Journal of Vegetation Science (2011).

- Pauchard, A. et al. Alien Plants Homogenise Protected Areas: Evidence
  from the Landscape and Regional Scales in South Central Chile. in
  Plant Invasions in Protected Areas (2013).

- Peyre, G. et al. VegPáramo, a flora and vegetation database for the
  Andean páramo. phytocoenologia (2015).

- Vibrans, A. C. et al. Insights from a large-scale inventory in the
  southern Brazilian Atlantic Forest. Scientia Agricola (2020).
