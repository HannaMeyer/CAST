test_that("trainDI works in default for a trained model", {
dat <- loaddata()
#...then calculate the DI of the trained model:
DI <- trainDI(model=dat$model)

#test threshold:
expect_equal(as.numeric(round(DI$threshold,5)), 0.38986)
# test summary statistics of the DI
expect_equal(as.numeric(colMeans(DI$train)),
             c(795.4426351,4.0277978,0.2577245))
})
