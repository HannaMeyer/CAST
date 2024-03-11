

test_that("AOA works in default: used with raster data and a trained model", {
  dat <- loaddata()
  # calculate the AOA of the trained model for the study area:
  AOA <- aoa(dat$studyArea, dat$model)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.38986)
  #test number of pixels within AOA:
  expect_equal(sum(values(AOA$AOA)==1,na.rm=TRUE), 2936)
  # test summary statistics of the DI
  expect_equal(as.vector(summary(values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1329  ", "Median :0.2052  ",
                 "Mean   :0.2858  ", "3rd Qu.:0.3815  ",
                "Max.   :4.4485  ", "NA's   :1993  "))
})


test_that("AOA works without a trained model", {
  dat <- loaddata()
  AOA <- aoa(dat$studyArea,train=dat$trainDat,variables=dat$variables)

  #test threshold:
  expect_equal(as.numeric(round(AOA$parameters$threshold,5)), 0.52872)
  #test number of pixels within AOA:
  expect_equal(sum(values(AOA$AOA)==1,na.rm=TRUE), 3377)
  # test summary statistics of the DI
  expect_equal(as.vector(summary(values(AOA$DI))),
               c("Min.   :0.0000  ", "1st Qu.:0.1759  ", "Median :0.2642  ",
                 "Mean   :0.3109  ", "3rd Qu.:0.4051  ",
                 "Max.   :2.6631  ", "NA's   :1993  "))
})

