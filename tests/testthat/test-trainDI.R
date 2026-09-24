source("data-fixture.R")

test_that("trainDI works in default for a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  #...then calculate the DI of the trained model:
  DI <- trainDI(model=dat$model, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(DI$threshold,5)), 0.38986)
  # test trainDI
  expect_equal(DI$trainDI, c(0.09043580, 0.14046341, 0.16584582, 0.57617177, 0.26840303,
                             0.14353894, 0.19768329, 0.24022059, 0.06832037, 0.29150668,
                             0.18471625, 0.57617177, 0.12344463, 0.09043580, 0.14353894,
                             0.26896008, 0.22713731, 0.24022059, 0.20388725, 0.06832037,
                             0.23604264, 0.20388725, 0.91513568, 0.09558666, 0.14046341,
                             0.16214832, 0.37107762, 0.16214832, 0.18471625, 0.12344463))
  # test summary statistics of the DI
  expect_equal(as.numeric(colMeans(DI$train)),
               c(795.4426351,4.0277978,0.2577245))
})

test_that("trainDI (with LPD = TRUE) works in default for a trained model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  #...then calculate the DI of the trained model:
  DI <- trainDI(model=dat$model, LPD = TRUE, verbose = F)

  #test threshold:
  expect_equal(as.numeric(round(DI$threshold,5)), 0.38986)
  #test trainLPD
  expect_identical(DI$trainLPD, as.integer(c(3, 4, 6, 0, 7,
                                             6, 2, 1, 5, 3,
                                             4, 0, 1, 2, 6,
                                             5, 4, 4, 5, 7,
                                             3, 4, 0, 2, 3,
                                             6, 1, 7, 3, 2)))
  # test summary statistics of the DI
  expect_equal(as.numeric(colMeans(DI$train)),
               c(795.4426351,4.0277978,0.2577245))
})

test_that("print and plot for trainDI run and return invisibly", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  DI <- trainDI(model = dat$model, verbose = FALSE)

  expect_no_error(print(DI))
  expect_invisible(print(DI))
  expect_s3_class(plot(DI), "ggplot")
})
