test_that("geodist works with points and polygon in geographic space", {

  data(splotdata)
  studyArea <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf")
  set.seed(1)
  folds <- data.frame("folds"=sample(1:3, nrow(splotdata), replace=TRUE))
  folds <- CreateSpacetimeFolds(folds, spacevar="folds", k=3)

  dist_geo <- geodist(x=splotdata,
                      modeldomain=studyArea,
                      cvfolds=folds$indexOut,
                      type = "geo")

  mean_sample2sample <- round(mean(dist_geo[dist_geo$what=="sample-to-sample","dist"]))
  mean_CV_distances <- round(mean(dist_geo[dist_geo$what=="CV-distances","dist"]))
  nrow_dist <- nrow(dist_geo)

  expect_equal(mean_sample2sample, 20321)
  expect_equal(mean_CV_distances, 25616)
  expect_equal(nrow_dist, 3410)


})

test_that("geodist works with points and polygon in feature space", {

  data(splotdata)
  studyArea <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf")
  set.seed(1)
  folds <- data.frame("folds"=sample(1:3, nrow(splotdata), replace=TRUE))
  folds <- CreateSpacetimeFolds(folds, spacevar="folds", k=3)
  predictors <- terra::rast(system.file("extdata","predictors_chile.tif", package="CAST"))

  dist_fspace <- geodist(x = splotdata,
                         modeldomain = predictors,
                         cvfolds=folds$indexOut,
                         type = "feature",
                         variables = c("bio_1","bio_12", "elev"))

  mean_sample2sample <- round(mean(dist_fspace[dist_fspace$what=="sample-to-sample","dist"]), 4)
  mean_CV_distances <- round(mean(dist_fspace[dist_fspace$what=="CV-distances","dist"]), 4)

  expect_equal(mean_sample2sample, 0.0843)
  expect_equal(mean_CV_distances, 0.1036)

})


test_that("geodist works space with points and preddata in geographic space", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")
  set.seed(1)
  ppoints <- suppressWarnings(sf::st_sample(aoi, 20, type="regular")) |>
    sf::st_set_crs("epsg:25832")

  set.seed(1)
  folds <- data.frame("folds"=sample(1:3, length(tpoints), replace=TRUE))
  folds <- CreateSpacetimeFolds(folds, spacevar="folds", k=3)

  dist_geo <- geodist(x=tpoints,
                      modeldomain=aoi,
                      preddata=ppoints,
                      type = "geo")

  mean_sample2sample <- round(mean(dist_geo[dist_geo$what=="sample-to-sample","dist"]), 4)
  mean_prediction_to_sample <- round(mean(dist_geo[dist_geo$what=="prediction-to-sample","dist"]), 4)

  expect_equal(mean_sample2sample, 1.4274)
  expect_equal(mean_prediction_to_sample, 2.9402)


})


test_that("geodist works with points and preddata in feature space", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")
  set.seed(1)
  ppoints <- suppressWarnings(sf::st_sample(aoi, 20, type="regular")) |>
    sf::st_set_crs("epsg:25832")

  raster <- terra::rast(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=10,
                        ymin=0, ymax=10, crs="epsg:25832", vals=1:100)

  dist <- geodist(x=tpoints,
                  modeldomain=raster,
                  preddata=ppoints,
                  type = "feature")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  mean_prediction_to_sample <- round(mean(dist[dist$what=="prediction-to-sample","dist"]), 4)

  expect_equal(mean_sample2sample, 0.4316)
  expect_equal(mean_prediction_to_sample, 0.8328)


})


test_that("geodist works with points and raster in geographic space", {

  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  raster <- terra::rast(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=10,
                        ymin=0, ymax=10, crs="epsg:25832", vals=1:100)

  dist <- geodist(x=tpoints,
                  modeldomain=raster,
                  type = "geo")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  expect_equal(mean_sample2sample, 1.4274)


})


test_that("geodist works with points and raster in feature space", {

  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  raster <- terra::rast(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=10,
                        ymin=0, ymax=10, crs="epsg:25832", vals=1:100)

  dist <- geodist(x=tpoints,
                  modeldomain=raster,
                  type = "feature")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  expect_equal(mean_sample2sample, 0.4316)


})


test_that("geodist works with points and stars raster in geographic space", {

  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  raster <- terra::rast(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=10,
                        ymin=0, ymax=10, crs="epsg:25832", vals=1:100) |>
    stars::st_as_stars()

  dist <- geodist(x=tpoints,
                  modeldomain=raster,
                  type = "feature")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  expect_equal(mean_sample2sample, 0.4316)


})



test_that("geodist works with points and test data in geographic space", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  set.seed(1)
  test_point <- suppressWarnings(sf::st_sample(aoi, 20, type="regular")) |>
    sf::st_set_crs("epsg:25832")

  dist <- geodist(x=tpoints,
                      modeldomain=aoi,
                      testdata=test_point,
                      type = "geo")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  mean_test_to_sample <- round(mean(dist[dist$what=="test-to-sample","dist"]), 4)

  expect_equal(mean_sample2sample, 1.4274)
  expect_equal(mean_test_to_sample, 2.9402)



})


test_that("geodist works with points and test data in feature space", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")

  raster <- terra::rast(nrows=10, ncols=10, nlyrs=1, xmin=0, xmax=10,
                        ymin=0, ymax=10, crs="epsg:25832", vals=1:100)

  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  set.seed(1)
  test_points <- suppressWarnings(sf::st_sample(aoi, 20, type="random")) |>
    sf::st_set_crs("epsg:25832")

  dist <- geodist(x=tpoints,
                  modeldomain=raster,
                  testdata = test_points,
                  type = "feature")

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  mean_test_to_sample <- round(mean(dist[dist$what=="test-to-sample","dist"]), 4)

  expect_equal(mean_sample2sample, 0.4316)
  expect_equal(mean_test_to_sample, 0.8783)


})


test_that("geodist works with categorical variables in feature space", {

  set.seed(1234)
  predictor_stack <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  predictors <- c("DEM","TWI", "NDRE.M", "Easting", "Northing", "fct")
  predictor_stack$fct <- factor(c(rep(LETTERS[1], terra::ncell(predictor_stack)/2),
                                  rep(LETTERS[2], terra::ncell(predictor_stack)/2)))

  predictor_stack <- predictor_stack[[predictors]]
  studyArea <- predictor_stack
  studyArea[!is.na(studyArea)] <- 1
  studyArea <- terra::as.polygons(studyArea, values = FALSE, na.all = TRUE) |>
    sf::st_as_sf() |>
    sf::st_union()

  pts <- clustered_sample(studyArea, 30, 5, 60)
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))
  pts <- terra::extract(predictor_stack, terra::vect(pts), ID=FALSE, bind=TRUE) |>
    sf::st_as_sf()

  test_pts <- clustered_sample(studyArea, 50, 5, 20)

  folds <- data.frame("folds"=sample(1:3, nrow(pts), replace=TRUE))
  folds <- CreateSpacetimeFolds(folds, spacevar="folds", k=3)

  sf::st_as_sf(terra::extract(predictor_stack,terra::vect(pts$geometry),bind=TRUE))

  dist <- geodist(x=pts,
                  modeldomain=predictor_stack,
                  type = "feature",
                  testdata = test_pts,
                  cvfolds = folds$indexOut)

  mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
  mean_prediction2sample <- round(mean(dist[dist$what=="prediction-to-sample","dist"]), 4)
  mean_test2sample <- round(mean(dist[dist$what=="test-to-sample","dist"]), 4)
  mean_CV_distance <- round(mean(dist[dist$what=="CV-distances","dist"]), 4)

  expect_equal(mean_sample2sample, 0.0459)
  expect_equal(mean_prediction2sample, 0.1625)
  expect_equal(mean_test2sample, 0.2358)
  expect_equal(mean_CV_distance, 0.0663)
})


test_that("geodist works in temporal space", {

dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
dat <- sf::st_as_sf(dat,coords=c("Easting","Northing"))
sf::st_crs(dat) <- 26911
trainDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2010,]
predictionDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2011,]
dist <- CAST::geodist(trainDat,preddata = predictionDat,type="time",time_unit="days")

mean_sample2sample <- round(mean(dist[dist$what=="sample-to-sample","dist"]), 4)
mean_prediction_to_sample <- round(mean(dist[dist$what=="prediction-to-sample","dist"]), 4)

expect_equal(mean_sample2sample, 0.02)
expect_equal(mean_prediction_to_sample, 194.7656)

dist <- CAST::geodist(trainDat,preddata = predictionDat,type="time",time_unit="hours")
mean_prediction_to_sample <- round(mean(dist[dist$what=="prediction-to-sample","dist"]), 2)
expect_equal(mean_prediction_to_sample, 4674.37)

})

test_that("geodist works in temporal space and with CV", {
  dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
  dat <- sf::st_as_sf(dat,coords=c("Easting","Northing"))
  sf::st_crs(dat) <- 26911
  trainDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2010,]
  predictionDat <- dat[dat$altitude==-0.3&lubridate::year(dat$Date)==2011,]
  trainDat$week <- lubridate::week(trainDat$Date)
  set.seed(100)
  cvfolds <- CreateSpacetimeFolds(trainDat,timevar = "week")

  dist <- CAST::geodist(trainDat,preddata = predictionDat,cvfolds = cvfolds$indexOut,
                        type="time",time_unit="days")

  mean_cv <- round(mean(dist[dist$what=="CV-distances","dist"]), 4)

    expect_equal(mean_cv,  2.4048)
}
)


