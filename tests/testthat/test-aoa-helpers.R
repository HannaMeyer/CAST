test_that(".validate_LPD accepts proportions and integers and rejects bad input", {
  expect_equal(CAST:::.validate_LPD(0.2, 100), 20)
  expect_equal(CAST:::.validate_LPD(1, 10), 10)
  expect_equal(CAST:::.validate_LPD(5, 10), 5L)

  expect_error(CAST:::.validate_LPD("a", 10), "maxLPD must be a number")
  expect_error(CAST:::.validate_LPD(0, 10), "maxLPD cannot be negative or equal to 0")
  expect_error(CAST:::.validate_LPD(0.001, 100), "proportion you provided for maxLPD is too small", fixed = FALSE)
  expect_error(CAST:::.validate_LPD(2.5, 10), "should be a whole number")
  expect_error(CAST:::.validate_LPD(11, 10), "cannot be bigger than the number of training samples")
})

test_that(".get_categorical_variables finds factor and character columns", {
  df <- data.frame(num = 1:3, f = factor(c("a","b","a")), c = c("x","y","x"), stringsAsFactors = FALSE)
  res <- CAST:::.get_categorical_variables(df, c("num","f","c"))
  expect_setequal(res, c("f","c"))
})

test_that(".drop_unknown_levels drops unused factor levels and maps unknowns to NA", {
  ref <- data.frame(cat = factor(c("a","b","a")), x = 1:3)
  qry <- data.frame(cat = c("a","c","b"), x = 4:6)
  out <- CAST:::.drop_unknown_levels(ref, qry, "cat")
  expect_true(is.factor(out$reference$cat))
  expect_equal(levels(out$reference$cat), c("a","b"))
  expect_true(is.factor(out$query$cat))
  expect_true(is.na(out$query$cat[2]))
})

test_that(".create_dummy_variables applies dummies correctly to query and weight", {
  skip_if_not_installed("caret")
  library(caret) # to avoid running into https://github.com/topepo/caret/issues/380

  ref <- data.frame(cat = factor(c("a", "b", "a")), v = 1:3, stringsAsFactors = FALSE)
  qry  <- data.frame(cat = c("a", "b", NA), v = 4:6, stringsAsFactors = FALSE)
  weight <- data.frame(v=1, cat = 5) # single-row weight for the categorical variable

  out <- CAST:::.create_dummy_variables(ref, qry, weight, "cat")

  # original categorical column removed
  expect_false("cat" %in% names(out$reference))
  expect_false("cat" %in% names(out$query))
  expect_false("cat" %in% names(out$weight))

  # dummy cols created and match the dummyVars columns
  dummy_cols <- setdiff(names(out$reference), names(ref))
  cols <- c("v", dummy_cols)
  expect_setequal(cols, names(out$reference))
  expect_setequal(cols, names(out$query))
  expect_setequal(cols, names(out$weight))

  # NA in query should have 0s in dummy columns
  expect_true(all(out$query[3, dummy_cols] == 0))

  # weight should have been expanded to the dummy columns
  expect_equal(as.numeric(out$weight[1, dummy_cols]), rep(as.numeric(weight$cat), length(dummy_cols)))
})

test_that(".convert_factors_to_dummy and .prepare_categorical_variables work together", {
  skip_if_not_installed("caret")
  library(caret) # to avoid running into https://github.com/topepo/caret/issues/380
  
  ref <- data.frame(a = factor(c("x","y")), b = 1:2)
  qry <- data.frame(a = c("x","z"), b = 3:4)
  res <- CAST:::.convert_factors_to_dummy(ref, qry, NULL, c("a"))
  expect_false("a" %in% names(res$reference))
  expect_false("a" %in% names(res$query))

  res2 <- CAST:::.prepare_categorical_variables(ref, qry, NULL, variables = c("a","b"))
  expect_true(is.list(res2))
  expect_true("catvars" %in% names(res2))
  expect_true("a" %in% res2$catvars)
})

test_that(".check_weights validates weights correctly", {
  vars <- c("x","y","z")
  w <- data.frame(x = 1, y = 2, z = 3)
  expect_equal(CAST:::.check_weights(w, vars), w)

  # list -> data.frame
  wlist <- as.list(w)
  expect_equal(as.numeric(CAST:::.check_weights(wlist, vars)[1, ]), as.numeric(w[1, ]))

  # missing variables -> error
  wbad <- data.frame(x = 1, y = 2)
  expect_error(CAST:::.check_weights(wbad, vars), "weights are not correctly specified")

  # multiple rows -> reset to ones with message
  wmulti <- data.frame(x = c(1,2), y = c(1,2), z = c(1,2))
  expect_message(wout <- CAST:::.check_weights(wmulti, vars), "weights are not correctly specified")
  expect_equal(as.numeric(wout[1, ]), rep(1, length(vars)))

  # negative weights -> set to zero with message
  wneg <- data.frame(x = -1, y = 0.5, z = 0.5)
  expect_message(wfixed <- CAST:::.check_weights(wneg, vars), "negative weights were set to 0")
  expect_true(all(wfixed >= 0))

  # weights sum to 0 -> error
  wzero <- data.frame(x = 0, y = 0, z = 0)
  expect_error(CAST:::.check_weights(wzero, vars), "weights sum to 0")
})

test_that(".prepare_weights respects useWeight flag and returns default weights", {
  vars <- c("a","b")
  w <- CAST:::.prepare_weights(NULL, NULL, vars, useWeight = FALSE)
  expect_true(is.data.frame(w))
  expect_equal(as.numeric(w[1, ]), rep(1, length(vars)))
  expect_equal(names(w), vars)
})

test_that(".apply_weights multiplies columns by weights and validates inputs", {
  df <- data.frame(a = 1:3, b = 2:4)
  w <- data.frame(a = 2, b = 3)
  out <- CAST:::.apply_weights(df, w)
  expect_equal(as.numeric(out$a), as.numeric(df$a * 2))
  expect_equal(as.numeric(out$b), as.numeric(df$b * 3))

  # NULL weight returns input
  expect_equal(CAST:::.apply_weights(df, NULL), df)

  # non-matching names -> error
  wbad <- data.frame(x = 1, y = 2)
  expect_error(CAST:::.apply_weights(df, wbad), "Weight columns do not match")

  # multi-row weight -> error
  wmulti <- data.frame(a = c(1,2), b = c(1,2))
  expect_error(CAST:::.apply_weights(df, wmulti), "Expected weight to be a single row")
})

test_that(".prepare_folds handles useCV flag, lists and vectors", {
  # useCV FALSE -> NULL with message
  expect_message(res0 <- CAST:::.prepare_folds(NULL, NULL, NULL, useCV = FALSE), "useCV is set to FALSE")
  expect_null(res0)

  # provide CVtrain as list and CVtest as vector
  CVtrain <- list(c(1,2), c(3,4))
  CVtest_vec <- c(2,1,2,1)
  res <- CAST:::.prepare_folds(NULL, CVtrain = CVtrain, CVtest = CVtest_vec, useCV = TRUE)
  expect_true(is.list(res))
  expect_true(all(c("CVtrain","CVtest") %in% names(res)))

  # mismatched lengths -> error
  CVtrain2 <- list(c(1), c(2), c(3))
  CVtest2 <- list(c(1), c(2))
  expect_error(CAST:::.prepare_folds(NULL, CVtrain2, CVtest2, useCV = TRUE), "should have the same number of folds")
})

test_that(".prepare_variables expands 'all' and validates variables exist in train", {
  train <- data.frame(a = 1, b = 2)
  vars <- CAST:::.prepare_variables("all", NULL, train)
  expect_setequal(vars, names(train))

  # specifying unknown variable errors
  expect_error(CAST:::.prepare_variables(c("a","c"), NULL, train), "some of the specified variables are not included in the training data")
})

test_that(".di_threshold returns a sensible robust threshold capped at max", {
  trainDI <- c(0.1, 0.2, 0.4, 1, 5)
  th <- CAST:::.di_threshold(trainDI)
  expected <- min(as.numeric(stats::quantile(trainDI, 0.75) + 1.5 * stats::IQR(trainDI)), max(trainDI))
  expect_equal(th, expected)
})

test_that(".prepare_weights accepts user-specified one-row data.frame", {
  vars <- c("a","b","c")
  user_w <- data.frame(a = 2, b = 0.5, c = 1)
  res <- CAST:::.prepare_weights(weight = user_w, variables = vars, useWeight = TRUE)
  expect_s3_class(res, "data.frame")
  expect_equal(as.numeric(res[1, ]), c(2, 0.5, 1))
  expect_equal(names(res), vars)
})