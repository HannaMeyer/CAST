
test_that("bss works with default arguments",{
  skip_on_cran()
  skip_on_os("mac", arch = "aarch64")
  skip_if_not_installed("randomForest")
  data("splotdata")
  splotdata = splotdata |> sf::st_drop_geometry()
  set.seed(1)
  selection = bss(predictors = splotdata[,6:12],
                  response = splotdata$Species_richness,
                  seed = 1,
                  verbose = FALSE,
                  ntree = 5,
                  tuneLength = 1)

  expect_identical(selection$selectedvars, c("bio_1", "bio_4", "bio_5", "bio_6", "bio_8", "bio_12"))
  expect_identical(selection$metric, "RMSE")
  expect_identical(selection$maximize, FALSE)
  expect_identical(round(selection$results$RMSE,2), 27.40)

})
