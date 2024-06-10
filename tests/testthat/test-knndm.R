test_that("kNNDM works with geographical coordinates and prediction points", {
  sf::sf_use_s2(TRUE)
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:4326")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:4326") |>
    sf::st_cast("POINT")
  set.seed(1)
  predpoints <- suppressWarnings(sf::st_sample(aoi, 20, type="regular")) |>
    sf::st_set_crs("epsg:4326")

  set.seed(1)
  kout <- knndm(tpoints, predpoints=predpoints, k=2, maxp=0.8)
  expect_identical(round(kout$W,1), 121095.2)
  expect_identical(kout$method, "hierarchical")
  expect_identical(kout$q, 3L)

})

test_that("kNNDM works with projected coordinates and prediction points", {
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")
  set.seed(1)
  predpoints <- sf::st_sample(aoi, 20, type="regular") |>
    sf::st_set_crs("epsg:25832")

  set.seed(1)
  kout <- knndm(tpoints, predpoints=predpoints, k=2, maxp=0.8, clustering = "kmeans")

  expect_identical(round(kout$W,4), 1.0919)
  expect_identical(kout$method, "kmeans")
  expect_identical(kout$q, 4L)

})

test_that("kNNDM works without crs and prediction points", {
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  tpoints <- sf::st_cast(sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))"), "POINT")
  set.seed(1)
  predpoints <- sf::st_sample(aoi, 20, type="regular")

  set.seed(1)
  kout <- suppressWarnings(knndm(tpoints, predpoints=predpoints, k=2, maxp=0.8))

  expect_identical(round(kout$W,6), 1.091896)
  expect_identical(kout$q, 3L)

  expect_warning(knndm(tpoints, predpoints=predpoints, k=2, maxp=0.8),
                 "Missing CRS in training or prediction points. Assuming projected CRS.")

})


test_that("kNNDM works with modeldomain and projected coordinates", {
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  set.seed(1)
  kout <- suppressMessages(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8, clustering = "kmeans"))

  expect_identical(round(kout$W,4), 1.2004)
  expect_identical(kout$method, "kmeans")
  expect_identical(kout$q, 4L)

  expect_message(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8, clustering = "kmeans"),
                 "1000 prediction points are sampled from the modeldomain")

})

test_that("kNNDM works with modeldomain and geographical coordinates", {
  sf::sf_use_s2(TRUE)
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:4326")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:4326") |>
    sf::st_cast("POINT")

  set.seed(1)
  kout <- suppressMessages(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8, clustering = "hierarchical"))

  expect_identical(round(kout$W,4), 133187.4275)
  expect_identical(kout$method, "hierarchical")
  expect_identical(kout$q, 3L)

  expect_message(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8, clustering = "hierarchical"),
                 "1000 prediction points are sampled from the modeldomain")

})

test_that("kNNDM works with modeldomain and no crs", {
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))") |>
    sf::st_cast("POINT")

  set.seed(1)
  kout <- suppressWarnings(suppressMessages(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8)))

  expect_identical(round(kout$W,4), 1.2004)
  expect_identical(kout$method, "hierarchical")
  expect_identical(kout$q, 3L)

  expect_message(suppressWarnings(knndm(tpoints, modeldomain = aoi, k=2, maxp=0.8)),
                 "1000 prediction points are sampled from the modeldomain")

})

test_that("kNNDM works when no clustering is present", {
  sf::sf_use_s2(TRUE)
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")

  set.seed(1)
  tpoints <- sf::st_sample(aoi, 10)

  set.seed(1)
  predpoints <- sf::st_sample(aoi, 20, type="regular")

  set.seed(1)
  kout <- suppressMessages(knndm(tpoints, predpoints = predpoints, k=2, maxp=0.8, clustering = "kmeans"))
  expect_equal(kout$q, "random CV")

  # for geographical coordinates
  set.seed(1)
  kout <- suppressMessages(knndm(sf::st_transform(tpoints,"epsg:4326"),
                                 predpoints = sf::st_transform(predpoints, "epsg:4326"),
                                 k=2, maxp=0.8, clustering = "hierarchical"))
  expect_equal(kout$q, "random CV")
})


test_that("kNNDM works with many points and different configurations", {
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  sample_area <- sf::st_as_sfc("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))", crs="epsg:25832")

  set.seed(1)
  tpoints <- sf::st_sample(sample_area, 100)

  set.seed(1)
  predpoints <- sf::st_sample(aoi, 1000)

  ks <- 2:10
  ps <- (1/ks)+0.1
  tune_grid <- data.frame(ks=ks, ps=ps)

  set.seed(1)
  kout <- apply(tune_grid, 1, function(j) {
    knndm(tpoints, predpoints=predpoints, k=j[[1]], maxp=j[[2]], clustering = "kmeans")
  })

  kout_W <- sapply(kout, function(x) round(x$W,3))
  kout_Gij <- sapply(kout, function(x) round(x$Gij[1],4))
  kout_Gjstar <- sapply(kout, function(x) round(x$Gjstar[1],4))

  w_expected <- c(2.184, 2.286, 2.468, 2.554, 2.570, 2.634, 2.678, 2.694, 2.688)
  Gij_expected <- rep(1.3886, length(w_expected))
  Gjstar_expected <- c(1.0981, 1.0981, 0.5400, 0.3812, 0.2505, 0.3812, 0.3099, 0.3099, 0.3812)

  expect_identical(round(kout_W,3), w_expected)
  expect_identical(round(kout_Gij,4), Gij_expected)
  expect_identical(round(Gjstar_expected,4), Gjstar_expected)

})


test_that("kNNDM recognizes erroneous input", {
  sf::sf_use_s2(TRUE)
  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  set.seed(1)
  predpoints <- sf::st_sample(aoi, 20)

  # maxp to small
  expect_error(knndm(tpoints, predpoints=predpoints, k=2, maxp=0.4))
  # k larger than number of tpoints
  expect_error(knndm(tpoints, predpoints=predpoints, k=20, maxp=0.8))
  # different crs of tpoints and predpoints
  expect_error(knndm(tpoints, predpoints=sf::st_transform(predpoints, "epsg:25833"), k=2, maxp=0.8))
  # different crs of tpoints and modeldomain
  expect_error(knndm(tpoints, modeldomain=sf::st_transform(aoi, "epsg:25833"), k=2, maxp=0.8))
  # using kmeans with geographical coordinates
  expect_error(knndm(sf::st_transform(tpoints,"epsg:4326"), predpoints=sf::st_transform(predpoints, "epsg:4326"),
                     clustering="kmeans"))
})

test_that("kNNDM yields the expected results with SpatRast modeldomain", {
  set.seed(1234)

  # prepare sample data
  data(cookfarm)
  dat <- terra::aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
                          by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- dat[,-1]
  pts <- sf::st_as_sf(pts,coords=c("Easting","Northing"))
  sf::st_crs(pts) <- 26911
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))

  knndm_folds <- knndm(pts, modeldomain = studyArea)
  expect_equal(as.numeric(knndm(pts, modeldomain = studyArea)$Gjstar[40]), 61.935505)
})


test_that("kNNDM works in feature space with kmeans clustering and raster as modeldomain", {
  set.seed(1234)

  # prepare sample data
  data(cookfarm)
  dat <- terra::aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
                          by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- dat[,-1]
  pts <- sf::st_as_sf(pts,coords=c("Easting","Northing"))
  sf::st_crs(pts) <- 26911
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))

  studyArea <- studyArea[[names(studyArea) %in% names(pts)]]
  train_points <- pts[,names(pts) %in% names(studyArea)]

  knndm_folds <- knndm(train_points, modeldomain = studyArea, space="feature", clustering = "kmeans")

  expect_equal(round(as.numeric(knndm_folds$Gjstar[40]),4), 0.2132)

})


test_that("kNNDM works in feature space with hierarchical clustering and raster as modeldomain", {
  set.seed(1234)

  # prepare sample data
  data(cookfarm)
  dat <- terra::aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
                          by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- dat[,-1]
  pts <- sf::st_as_sf(pts,coords=c("Easting","Northing"))
  sf::st_crs(pts) <- 26911
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))

  studyArea <- studyArea[[names(studyArea) %in% names(pts)]]
  tpoints <- pts[,names(pts) %in% names(studyArea)]

  knndm_folds <- knndm(tpoints, modeldomain = studyArea, space="feature", clustering = "hierarchical")

  expect_equal(round(as.numeric(knndm_folds$Gjstar[40]),4), 0.2132)


})

test_that("kNNDM works in feature space with clustered training points", {
  set.seed(1234)

  data(splotdata)
  splotdata <- splotdata[splotdata$Country == "Chile",]

  predictors <- c("bio_1", "bio_4", "bio_5", "bio_6",
                  "bio_8", "bio_9", "bio_12", "bio_13",
                  "bio_14", "bio_15", "elev")
  trainDat <- sf::st_drop_geometry(splotdata)
  predictors_sp <- terra::rast(system.file("extdata", "predictors_chile.tif",package="CAST"))

  knndm_folds <- knndm(trainDat[,predictors], modeldomain = predictors_sp, space = "feature",
                       clustering="kmeans", k=4, maxp=0.8)


  expect_equal(round(as.numeric(knndm_folds$Gjstar[40]),4), 0.8287)

})


test_that("kNNDM works in feature space with categorical variables and predpoints", {
  set.seed(1234)

  # prepare sample data
  data(cookfarm)
  dat <- terra::aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
                          by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- dat[,-1]
  pts <- sf::st_as_sf(pts,coords=c("Easting","Northing"))
  sf::st_crs(pts) <- 26911
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))

  studyArea <- studyArea[[names(studyArea) %in% names(pts)]]

  prediction_points <- terra::spatSample(studyArea, 1000, "regular")
  train_points <- pts[,names(pts) %in% names(studyArea)]

  prediction_points$fct <- factor(sample(LETTERS[1:4], nrow(prediction_points), replace=TRUE))
  train_points$fct <- factor(sample(LETTERS[1:4], nrow(pts), replace=TRUE))


  knndm_folds <- knndm(tpoints=train_points, predpoints = prediction_points,
                       space="feature", clustering = "hierarchical")

  expect_equal(round(as.numeric(knndm_folds$Gjstar[40]),3), 0.057)

})


test_that("kNNDM works in feature space with clustered training points, categorical features ", {
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
  pts <- terra::extract(predictor_stack, terra::vect(pts), ID=FALSE)

  knndm_folds_kproto <- knndm(tpoints=pts, modeldomain = predictor_stack, space="feature", clustering = "kmeans")
  knndm_folds_hclust <- knndm(tpoints=pts, modeldomain = predictor_stack, space="feature", clustering = "hierarchical")

  expect_equal(round(as.numeric(knndm_folds_kproto$Gjstar[20]),3), 0.077)
  expect_equal(round(as.numeric(knndm_folds_hclust$Gjstar[20]),3), 0.078)

})


test_that("kNNDM works in feature space with Mahalanobis distance", {
  data(splotdata)
  splotdata <- splotdata[splotdata$Country == "Chile",]

  predictors <- c("bio_1", "bio_4", "bio_5", "bio_6",
                  "bio_8", "bio_9", "bio_12", "bio_13",
                  "bio_14", "bio_15", "elev")
  trainDat <- sf::st_drop_geometry(splotdata)
  predictors_sp <- terra::rast(system.file("extdata", "predictors_chile.tif",package="CAST"))

  set.seed(1234)
  knndm_folds <- knndm(trainDat[,predictors], modeldomain = predictors_sp, space = "feature",
                       clustering="kmeans", k=4, maxp=0.8, useMD=TRUE)

  expect_equal(round(as.numeric(knndm_folds$Gjstar[40]),4), 1.1258)

})


test_that("kNNDM works in feature space with Mahalanobis distance without clustering", {
  set.seed(1234)

  # prepare sample data
  data(cookfarm)
  dat <- terra::aggregate(cookfarm[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
                          by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- dat[,-1]
  pts <- sf::st_as_sf(pts,coords=c("Easting","Northing"))
  sf::st_crs(pts) <- 26911
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
  pts <- sf::st_transform(pts, crs = sf::st_crs(studyArea))

  studyArea <- studyArea[[names(studyArea) %in% names(pts)]]
  train_points <- pts[,names(pts) %in% names(studyArea)]

  expect_message(knndm(train_points, modeldomain = studyArea, space="feature", clustering = "kmeans", useMD = TRUE),
                 "Gij <= Gj; a random CV assignment is returned")

  expect_message(knndm(train_points, modeldomain = studyArea, space="feature", clustering = "hierarchical", useMD = TRUE),
                 "Gij <= Gj; a random CV assignment is returned")

})
