# Cookfarm soil logger data

spatio-temporal data of soil properties and associated predictors for
the Cookfarm in Washington, USA. The data are a subset of the cookfarm
dataset provided with the [GSIF
package](https://CRAN.R-project.org/package=GSIF).

## Usage

``` r
data(cookfarm)
```

## Format

A sf data.frame with 128545 rows and 17 columns:

- SOURCEID:

  ID of the logger

- VW:

  Response Variable - Soil Moisture

- altitude:

  Measurement depth of VW

- Date, cdata:

  Measurement Date, Cumulative Date

- Easting, Northing:

  Location Coordinates (EPSG:26911)

- DEM, TWI, NDRE.M, NDRE.Sd, Precip_wrcc, MaxT_wrcc, MinT_wrcc,
  Precip_cum:

  Predictor Variables

## References

- Gash et al. 2015 - Spatio-temporal interpolation of soil water,
  temperature, and electrical conductivity in 3D + T: The Cook Agronomy
  Farm data set
  [doi:10.1016/j.spasta.2015.04.001](https://doi.org/10.1016/j.spasta.2015.04.001)

- Meyer et al. 2018 - Improving performance of spatio-temporal machine
  learning models using forward feature selection and target-oriented
  validation
  [doi:10.1016/j.envsoft.2017.12.001](https://doi.org/10.1016/j.envsoft.2017.12.001)
