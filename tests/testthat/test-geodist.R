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
  # can't be tested for prediction-to-sample, which are sampled slightly different in each run
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
  # can't be tested for prediction-to-sample, which are sampled slightly different in each run

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
