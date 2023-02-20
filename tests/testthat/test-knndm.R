test_that("kNNDM works with geographical coordinates and prediction points", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:4326")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:4326") |>
    sf::st_cast("POINT")
  set.seed(1)
  ppoints <- suppressWarnings(sf::st_sample(aoi, 20, type="regular")) |>
    sf::st_set_crs("epsg:4326")

  set.seed(1)
  kout <- knndm(tpoints, ppoints=ppoints, k=2, maxp=0.8)
  expect_identical(round(kout$W,1), 121095.2)
  expect_identical(kout$method, "hierarchical")
  expect_identical(kout$q, 3L)

})

test_that("kNNDM works with projected coordinates and prediction points", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")
  set.seed(1)
  ppoints <- sf::st_sample(aoi, 20, type="regular") |>
    sf::st_set_crs("epsg:25832")

  set.seed(1)
  kout <- knndm(tpoints, ppoints=ppoints, k=2, maxp=0.8, clustering = "kmeans")

  expect_identical(round(kout$W,4), 1.0919)
  expect_identical(kout$method, "kmeans")
  expect_identical(kout$q, 4L)

})

test_that("kNNDM works without crs and prediction points", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  tpoints <- sf::st_cast(sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))"), "POINT")
  set.seed(1)
  ppoints <- sf::st_sample(aoi, 20, type="regular")

  set.seed(1)
  kout <- suppressWarnings(knndm(tpoints, ppoints=ppoints, k=2, maxp=0.8))

  expect_identical(round(kout$W,6), 1.091896)
  expect_identical(kout$q, 3L)

  expect_warning(knndm(tpoints, ppoints=ppoints, k=2, maxp=0.8),
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

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")

  set.seed(1)
  tpoints <- sf::st_sample(aoi, 10)

  set.seed(1)
  ppoints <- sf::st_sample(aoi, 20, type="regular")

  set.seed(1)
  kout <- suppressMessages(knndm(tpoints, ppoints = ppoints, k=2, maxp=0.8, clustering = "kmeans"))
  expect_equal(kout$q, "random CV")

  # for geographical coordinates
  set.seed(1)
  kout <- suppressMessages(knndm(sf::st_transform(tpoints,"epsg:4326"),
                                 ppoints = sf::st_transform(ppoints, "epsg:4326"),
                                 k=2, maxp=0.8, clustering = "hierarchical"))
  expect_equal(kout$q, "random CV")
})


test_that("kNNDM works with many points and different configurations", {

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  sample_area <- sf::st_as_sfc("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))", crs="epsg:25832")

  set.seed(1)
  tpoints <- sf::st_sample(sample_area, 100)

  set.seed(1)
  ppoints <- sf::st_sample(aoi, 1000)

  ks <- 2:10
  ps <- (1/ks)+0.1
  tune_grid <- data.frame(ks=ks, ps=ps)

  set.seed(1)
  kout <- apply(tune_grid, 1, function(j) {
    knndm(tpoints, ppoints=ppoints, k=j[[1]], maxp=j[[2]], clustering = "kmeans")
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

  aoi <- sf::st_as_sfc("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="epsg:25832")
  tpoints <- sf::st_as_sfc("MULTIPOINT ((1 1), (1 2), (2 2), (2 3), (1 4), (5 4))", crs="epsg:25832") |>
    sf::st_cast("POINT")

  set.seed(1)
  ppoints <- sf::st_sample(aoi, 20)

  # maxp to small
  expect_error(knndm(tpoints, ppoints=ppoints, k=2, maxp=0.4))
  # k larger than number of tpoints
  expect_error(knndm(tpoints, ppoints=ppoints, k=20, maxp=0.8))
  # different crs of tpoints and ppoints
  expect_error(knndm(tpoints, ppoints=sf::st_transform(ppoints, "epsg:25833"), k=2, maxp=0.8))
  # different crs of tpoints and modeldomain
  expect_error(knndm(tpoints, modeldomain=sf::st_transform(aoi, "epsg:25833"), k=2, maxp=0.8))
  # using kmeans with geographical coordinates
  expect_error(knndm(sf::st_transform(tpoints,"epsg:4326"), ppoints=sf::st_transform(ppoints, "epsg:4326"),
                     clustering="kmeans"))
})

