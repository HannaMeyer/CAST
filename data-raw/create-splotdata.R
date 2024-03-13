## This script creates the example dataset "splotdata" of the CAST package.
## It downloads splotopen data points and associated worldclim predictors for South America.
## A lower resolution predictor stack (terra rast) is also created for Chile.
## For more information, please check out the Book Chapter and Repository CAST4Ecology

library(geodata)
library(rnaturalearth)
library(terra)
library(sf)
library(tidyverse)
library(geodata)


##### Download Predictors --------------------------------
## Warning: This downloads ~ 1 GB of data
dir.create("data-raw/raw/")

wcf = geodata::worldclim_global(var = "bio", path = "data-raw/raw/", res = 0.5)
wc = geodata::worldclim_global(var = "bio", path = "data-raw/raw/", res = 5)
elevf = geodata::elevation_global(res = 0.5, path = "data-raw/raw/")
elev = geodata::elevation_global(res = 5, path = "data-raw/raw/")

wcf = c(wcf, elevf)
wc = c(wc, elev)

##### Download sPlotOpen -------------------------------------
if(!file.exists("data-raw/raw/splotopen")){
  download.file("https://idata.idiv.de/ddm/Data/DownloadZip/3474?version=5779", destfile = "data-raw/raw/splotopen.zip")
  unzip(zipfile = "data-raw/raw/splotopen.zip", exdir = "data-raw/raw/splotopen")
  unzip(zipfile = "data-raw/raw/splotopen/sPlotOpen.RData(2).zip", exdir = "data-raw/raw/splotopen")
}



##### Clean up and save necessary files ----------------------------------
# define region: all of south america
region = rnaturalearth::ne_countries(continent = "South America", returnclass = "sf", scale = 110)


# Predictor clean up
wc = crop(wc, region)
names(wc) = names(wc) |> str_remove(pattern = "wc2.1_5m_")
p = c("bio_1", "bio_4", "bio_5", "bio_6", "bio_8", "bio_9", "bio_12", "bio_13", "bio_14", "bio_15", "elev")
wc = wc[[p]]

# worldclim in full resolution for extracting the training data
wcf = crop(wcf, region)
names(wcf) = names(wcf) |> str_remove(pattern = "wc2.1_30s_")
wcf = wcf[[p]]
wcf$lat = terra::init(wcf, "y")
wcf$lon = terra::init(wcf, "x")


# Gather Response Variable: sPlotOpen Species Richness for South America
## see Appendix 1 of https://doi.org/10.1111/geb.13346
load("data-raw/raw/splotopen/sPlotOpen.RData")

splot = header.oa |>
  #filter(Resample_1 == TRUE) |>
  filter(Continent == "South America") |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  left_join(CWM_CWV.oa |> select(c("PlotObservationID", "Species_richness"))) |>
  select(c("PlotObservationID", "GIVD_ID", "Country", "Biome",
           "Species_richness")) |>
  na.omit()

# extract predictor values and attach to response
splot = terra::extract(wcf, splot, ID = FALSE, bind = TRUE) |>
  st_as_sf() |>
  na.omit()


# only keep unique locations
## some reference sample locations are in the same predictor stack pixel
## this can lead to erroneous models and misleading validations
splotdata = splot[!duplicated(c(splot$lat, splot$lon)),]
splotdata = splotdata |> na.omit()
splotdata$lat = NULL
splotdata$lon = NULL


# save splotdata
splotdata$Biome = droplevels(splotdata$Biome)
save(splotdata, file = "data/splotdata.rda", compress = "xz")

## save predictors for chile
chile = rnaturalearth::ne_countries(country = "Chile", returnclass = "sf")
wc = crop(wc, chile)
writeRaster(wc, "inst/extdata/predictors_chile.tif", datatype = "INT2S", overwrite = TRUE)


## Remove downloaded data
unlink("data-raw/raw", recursive = TRUE)
