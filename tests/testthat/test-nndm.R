test_that("Valid range of phi", {
  set.seed(1234)
  poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                     byrow=TRUE)))
  poly_sfc <- sf::st_sfc(poly)
  tpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "random")
  predpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "regular")

  expect_error(nndm(tpoints_sfc, predpoints = predpoints_sfc, phi = -1),
               "phi must be positive or set to 'max'.")
})

test_that("NNDM detects wrong data and geometry types", {
  set.seed(1234)
  poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                     byrow=TRUE)))
  poly_sfc <- sf::st_sfc(poly)
  tpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "random")
  predpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "regular")

  # tpoints
  expect_error(suppressWarnings(nndm(1, predpoints = predpoints_sfc)),
               "tpoints must be a sf/sfc object.")
  expect_error(nndm(poly, predpoints = predpoints_sfc),
               "tpoints must be a sf/sfc object.")
  expect_error(nndm(sf::st_sfc(poly), predpoints = predpoints_sfc),
               "tpoints must be a sf/sfc point object.")
  # predpoints
  expect_error(suppressWarnings(nndm(tpoints_sfc, predpoints = 1)),
               "predpoints must be a sf/sfc object.")
  expect_error(nndm(tpoints_sfc, predpoints = poly),
               "predpoints must be a sf/sfc object.")
  expect_error(nndm(tpoints_sfc, predpoints = poly_sfc),
               "predpoints must be a sf/sfc point object.")

  # model domain
  expect_error(suppressWarnings(nndm(tpoints_sfc, modeldomain = 1)),
               "modeldomain must be a sf/sfc object or a 'SpatRaster' object.")
  expect_error(nndm(tpoints_sfc, modeldomain = predpoints_sfc),
               "modeldomain must be a sf/sfc polygon object.")
})

test_that("NNDM detects different CRS in inputs", {
  sf::sf_use_s2(TRUE)
  set.seed(1234)
  poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                     byrow=TRUE)))
  poly_sfc <- sf::st_sfc(poly)
  tpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "random")
  predpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "regular")

  tpoints_sfc_4326 <- sf::st_set_crs(tpoints_sfc, 4326)
  tpoints_sfc_3857 <- sf::st_set_crs(tpoints_sfc, 3857)
  predpoints_sfc_4326 <- sf::st_set_crs(predpoints_sfc, 4326)
  predpoints_sfc_3857 <- sf::st_set_crs(predpoints_sfc, 3857)
  poly_sfc_4326 <- sf::st_set_crs(poly_sfc, 4326)

  # tests
  expect_error(nndm(tpoints_sfc_3857, predpoints = predpoints_sfc),
               "tpoints and predpoints must have the same CRS")
  expect_error(nndm(tpoints_sfc_3857, modeldomain = poly_sfc_4326),
               "tpoints and modeldomain must have the same CRS")
})



test_that("NNDM yields the expected results for all data types", {
  set.seed(1234)
  poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                     byrow=TRUE)))
  poly_sfc <- sf::st_sfc(poly)
  tpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "random")
  predpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "regular")

  # tpoints, predpoints
  expect_equal(as.numeric(nndm(tpoints_sfc, predpoints = tpoints_sfc)$Gjstar[1]), 3.7265881)
  # tpoints, modeldomain
  expect_equal(as.numeric(nndm(tpoints_sfc, modeldomain = poly_sfc)$Gjstar[5]), 4.9417614)
  # change phi
  expect_equal(as.numeric(nndm(tpoints_sfc, predpoints = tpoints_sfc, phi = 10)$Gjstar[10]), 4.8651321)
  # change min_train
  expect_equal(as.numeric(nndm(tpoints_sfc, predpoints = tpoints_sfc, phi = 20, min_train = 0.2)$Gjstar[15]), 3.466861)
  # length checks
  expect_equal(length(nndm(tpoints_sfc, predpoints = tpoints_sfc)$Gjstar), length(tpoints_sfc))
  expect_equal(length(nndm(tpoints_sfc, predpoints = tpoints_sfc)$Gi), length(tpoints_sfc))
  expect_gt(length(nndm(tpoints_sfc, modeldomain = poly_sfc)$Gij), 900)
})

test_that("NNDM yields the expected results for all CRS", {
  sf::sf_use_s2(TRUE)
  set.seed(1234)
  poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                     byrow=TRUE)))
  poly_sfc <- sf::st_sfc(poly)
  tpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "random")
  predpoints_sfc <- sf::st_sample(poly_sfc, 50, type = "regular")

  # Projected
  tpoints_3857 <- sf::st_set_crs(tpoints_sfc, 3857)
  predpoints_3857 <- sf::st_set_crs(predpoints_sfc, 3857)
  expect_equal(as.numeric(nndm(tpoints_3857, predpoints = predpoints_3857, phi = 10)$Gjstar[20]), 3.2921498)

  # Geographic
  tpoints_sf_4326 <- sf::st_set_crs(tpoints_sfc, 4326)
  predpoints_sf_4326 <- sf::st_set_crs(predpoints_sfc, 4326)
  expect_equal(as.numeric(nndm(tpoints_sf_4326, predpoints = predpoints_sf_4326, phi = 1000000)$Gjstar[20]), 355614.94)
})

test_that("NNDM yields the expected results with SpatRast modeldomain", {
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

  nndm_folds <- nndm(pts, modeldomain = studyArea, phi = 150)
  expect_equal(as.numeric(nndm(pts, modeldomain = studyArea, phi = 150)$Gjstar[5]), 63.828663)

})
