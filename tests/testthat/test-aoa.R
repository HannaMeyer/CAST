source("data-fixture.R")

test_that("AOA works in default: used with raster data and a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 2936)
  # test trainDI
  expect_equal(AOA$parameters$trainDI, c(0.09043580, 0.14046341, 0.16584582, 0.57617177, 0.26840303,
                                         0.14353894, 0.19768329, 0.24022059, 0.06832037, 0.29150668,
                                         0.18471625, 0.57617177, 0.12344463, 0.09043580, 0.14353894,
                                         0.26896008, 0.22713731, 0.24022059, 0.20388725, 0.06832037,
                                         0.23604264, 0.20388725, 0.91513568, 0.09558666, 0.14046341,
                                         0.16214832, 0.37107762, 0.16214832, 0.18471625, 0.12344463))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI)))[1:6],
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                "Max.   :4.4485  "))
})

test_that("AOA works with a stars object", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("stars")
  dat <- loaddata()
  studyArea_stars <- stars::st_as_stars(dat$studyArea)
  AOA <- aoa(studyArea_stars, dat$model, LPD = TRUE, verbose = F)
  expect_true(inherits(AOA, "aoa"))
  expect_true(inherits(AOA$DI, "stars"))
  expect_true(inherits(AOA$AOA, "stars"))
  expect_true(inherits(AOA$LPD, "stars"))
})


test_that("AOA works without a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI)))[1:6],
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  "))
})

test_that("AOA (including LPD) works with raster data and a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 2936)
  #test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(3, 4, 6, 0, 7,
                                          6, 2, 1, 5, 3,
                                          4, 0, 1, 2, 6,
                                          5, 4, 4, 5, 7,
                                          3, 4, 0, 2, 3,
                                          6, 1, 7, 3, 2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI)))[1:6],
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                 "Max.   :4.4485  "))
})


test_that("AOA (inluding LPD) works without a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables, LPD = TRUE, maxLPD = 1, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(terra::values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test trainLPD
  expect_equal(AOA$parameters$trainLPD, c(7, 9, 12, 1, 12,
                                          12, 4, 2, 8, 10,
                                          6, 1, 3,4, 11,
                                          9, 9, 7, 5, 5,
                                          6, 5, 0, 5, 9,
                                          8, 4, 11, 3,2))
  # test summary statistics of the DI
  expect_equal(as.vector(summary(terra::values(AOA$DI)))[1:6],
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  "))
})


test_that("AOA raises warnings with deprecated parameters in verbose mode", {
  skip_on_cran()
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  expect_warning(
    AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = T, parallel = TRUE, cores = 2),
    "The 'parallel' and 'cores' parameters are deprecated.")
  expect_true(inherits(AOA, "aoa"))
  expect_warning(
    AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = T, method = "euclidean", algorithm = "brute"),
    "The 'method' and 'algorithm' parameters are deprecated.")
  expect_true(inherits(AOA, "aoa"))
})



test_that("print and plot for aoa run and return invisibly", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea, dat$model, verbose = FALSE)

  expect_no_error(print(AOA))
  expect_invisible(print(AOA))
  expect_s3_class(plot(AOA, samplesize = 10), "ggplot")
})

test_that("errorProfiles works for aoa objects (DI)", {
  skip_on_cran()
  skip_on_os("mac", arch = "aarch64")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("scam")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea, dat$model, verbose = FALSE)

  err_model <- errorProfiles(dat$model, AOA, variable = "DI")
  expect_s3_class(err_model, "errorModel")
  expect_true(is.numeric(attr(err_model, "AOA_threshold")))
})

test_that("errorProfiles works for aoa objects (LPD)", {
  skip_on_cran()
  skip_on_os("mac", arch = "aarch64")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("scam")
  dat <- loaddata()
  AOA <- aoa(dat$studyArea, dat$model, LPD = TRUE, maxLPD = 1, verbose = FALSE)

  err_model <- errorProfiles(dat$model, AOA, variable = "LPD")
  expect_s3_class(err_model, "errorModel")
  expect_true(is.numeric(attr(err_model, "AOA_threshold")))
})

test_that("AOA masks unknown factor levels in newdata (catvars) as NA", {
  # create simple training data with a factor predictor
  set.seed(42)
  train <- data.frame(x = rnorm(20), y = rnorm(20))
  train$fac <- factor(sample(c("a", "b"), nrow(train), replace = TRUE))

  # newdata contains an unseen level "c" which should be masked as NA
  newdata <- data.frame(x = c(0, 1, -1), y = c(0, 1, -1), fac = factor(c("a", "c", "b")))

  AOA <- aoa(newdata, train = train, variables = c("x", "y", "fac"), verbose = FALSE)
  
  expect_s3_class(AOA, "aoa")
  expect_equal(length(AOA$DI), nrow(newdata))
  expect_all_true(AOA$DI > 0) # DI should be > 0 for all rows
})

test_that("AOA works when newdata is a SpatRaster with unseen factor levels", {
  # training data with factor levels a,b
  train <- data.frame(x = c(0, 1, -1), y = c(0, 1, -1))
  train$fac <- factor(c("a", "b", "a"))

  # create a small SpatRaster with three cells and three layers: x, y, fac
  r_x <- terra::rast(nrows = 1, ncols = 3)
  terra::values(r_x) <- c(0, 1, -1)
  names(r_x) <- "x"

  r_y <- terra::rast(nrows = 1, ncols = 3)
  terra::values(r_y) <- c(0, 1, -1)
  names(r_y) <- "y"

  fac <- terra::rast(nrows = 1, ncols = 3)
  # values 1,2,3 correspond to levels a,b,c where c is unseen in train
  terra::values(fac) <- c(1, 3, 2)
  fac <- terra::as.factor(fac)
  levels(fac) <- data.frame(ID = 1:3, fac = c("a", "b", "c"))
  names(fac) <- "fac"

  new_r <- c(r_x, r_y, fac)

  AOA <- aoa(new_r, train = train, variables = c("x", "y", "fac"), verbose = FALSE)

  expect_s3_class(AOA, "aoa")
  expect_true(inherits(AOA$DI, "SpatRaster"))
  vals <- terra::values(AOA$DI)
  expect_equal(length(vals), terra::ncell(new_r))
  # ensure DI was computed (no NA / NaN values)
  expect_all_true(is.finite(vals[ ,1]))
})

