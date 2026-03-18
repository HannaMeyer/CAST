loaddata <- function() {
  # prepare sample data:
  data("cookfarm")
  dat <- aggregate(cookfarm[,c("VW","Easting","Northing")],by=list(as.character(cookfarm$SOURCEID)),mean)
  pts <- sf::st_as_sf(dat,coords=c("Easting","Northing"))
  pts$ID <- 1:nrow(pts)
  set.seed(100)
  pts <- pts[1:30,]
  studyArea <- terra::rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))[[1:8]]
  trainDat <- terra::extract(studyArea,pts,na.rm=FALSE)
  trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")

  # train a model:
  set.seed(100)
  variables <- c("DEM","NDRE.Sd","TWI")
  model <- caret::train(trainDat[,which(names(trainDat)%in%variables)],
                        trainDat$VW, method="rf", importance=TRUE, tuneLength=1,
                        trControl=caret::trainControl(method="cv",number=5,savePredictions=T))


  data <- list(
    studyArea = studyArea,
    trainDat = trainDat,
    variables = variables,
    model = model
  )

  return(data)
}

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

test_that("aoa_get_train returns model training data", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  tr <- aoa_get_train(dat$model)
  expect_equal(as.data.frame(dat$model$trainingData), as.data.frame(tr))
})

test_that("aoa_get_weights returns sensible weights for a caret model", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  w <- aoa_get_weights(dat$model, dat$variables)
  expect_s3_class(w, "data.frame")
  expect_true(all(dat$variables %in% names(w)))
  expect_true(sum(unlist(w)) > 0)
})

test_that("aoa_get_weights falls back to default when varImp extraction fails", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()

  # corrupt the trained model so varImp/importance will fail
  bad_model <- dat$model
  bad_model$finalModel <- NULL

  w <- aoa_get_weights(bad_model, dat$variables)

  expect_s3_class(w, "data.frame")
  expect_true(all(dat$variables %in% names(w)))
  expect_equal(as.numeric(w[1, ]), rep(1, length(dat$variables)))
})

test_that("user_weights accepts correct input and falls back on defaults for bad input", {
  vars <- c("a","b")
  good <- data.frame(a = 0.2, b = 0.8)
  res_good <- user_weights(good, vars)
  expect_equal(as.numeric(res_good[1,]), c(0.2, 0.8))

  bad <- data.frame(x = 1, y = 2)
  expect_message(res_bad <- user_weights(bad, vars), "variable weights are not correctly specified")
  expect_equal(as.numeric(res_bad[1,]), c(1,1))

  neg <- data.frame(a = -1, b = 0.5)
  expect_message(res_neg <- user_weights(neg, vars), "negative weights were set to 0")
  expect_equal(as.numeric(res_neg[1,]), c(0, 0.5))
})

test_that("aoa_categorial_train expands factor variables and adjusts weights", {
  train <- data.frame(Cat = factor(c("a","b","a")), Num = c(1,2,3))
  vars <- c("Cat","Num")
  weight <- data.frame(Cat = 0.5, Num = 1)
  expect_message(res <- aoa_categorial_train(train, vars, weight), "warning: predictors contain categorical variables")
  # categorical var reported
  expect_equal(res$catvars, "Cat")
  # original categorical column removed
  expect_false("Cat" %in% names(res$train))
  # numeric column preserved and dummy columns added (ncol = original -1 + nlevels)
  expect_equal(ncol(res$train), ncol(train) - 1 + length(levels(train$Cat)))
  # weight names match train column names
  expect_equal(names(res$weight), names(res$train))

  # when no categorical variables, returns unchanged structure
  train2 <- data.frame(Num = c(1,2,3))
  vars2 <- c("Num")
  weight2 <- data.frame(Num = 1)
  res2 <- aoa_categorial_train(train2, vars2, weight2)
  expect_equal(res2$catvars, character(0))
  expect_equal(res2$train, train2)
  expect_equal(res2$weight, weight2)
})

test_that("aoa_get_folds extracts folds from caret model and respects useCV flag", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  # when using model and default useCV = TRUE, folds are extracted from model
  folds <- aoa_get_folds(dat$model, CVtrain = NULL, CVtest = NULL, useCV = TRUE)
  expect_true(is.list(folds))
  expect_true(length(folds) == 2)
  # CVtrain and CVtest should match model control index structures
  expect_identical(folds[[1]], dat$model$control$index)
  expect_identical(folds[[2]], dat$model$control$indexOut)

  # when useCV is FALSE, folds should be NULL and a message is emitted
  expect_message(folds2 <- aoa_get_folds(dat$model, CVtrain = NULL, CVtest = NULL, useCV = FALSE),
                 "useCV is set to FALSE")
  expect_equal(folds2[[1]], NULL)
  expect_equal(folds2[[2]], NULL)
})


test_that("aoa_get_folds restructures CVtest vector into lists and derives complements", {
  # when model is NA and CVtest is a vector of fold ids it should be restructured
  CVtest_vec <- c(1,1,2,2,3,3)
  folds <- aoa_get_folds(NULL, CVtrain = NULL, CVtest = CVtest_vec, useCV = TRUE)
  expect_type(folds, "list")
  CVtest_list <- folds[[2]]
  CVtrain_list <- folds[[1]]
  expect_true(is.list(CVtest_list))
  expect_true(is.list(CVtrain_list))
  expect_equal(CVtest_list[[1]], which(CVtest_vec == 1))
  # complements: none of the indices in the train list element should be in the corresponding test element
  expect_true(all(!CVtrain_list[[1]] %in% CVtest_list[[1]]))
})


test_that("aoa_get_variables returns all model variables when 'all' is requested", {
  skip_if_not_installed("randomForest")
  dat <- loaddata()
  vars <- aoa_get_variables("all", dat$model, dat$trainDat)
  # should equal the variables we passed when training (order may differ)
  expect_equal(sort(vars), sort(dat$variables))
})

