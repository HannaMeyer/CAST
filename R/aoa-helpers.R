################################################################################
# Input validation
################################################################################
.validate_LPD <- function(maxLPD, n_samples) {
  if (!inherits(maxLPD, "numeric") && !inherits(maxLPD, "integer")) {
    stop("maxLPD must be a number. Either define a number between 0 and 1 to use a proportion of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
  }
  if (maxLPD <= 0) {
    stop("maxLPD cannot be negative or equal to 0. Either define a number between 0 and 1 to use a proportion of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
  }

  is_proportion <- maxLPD <= 1

  if (is_proportion) {
    maxLPD <- round(maxLPD * n_samples)
    if (maxLPD <= 1) {
      stop("The proportion you provided for maxLPD is too small.")
    }
  } else {
    if (maxLPD %% 1 != 0) {
      stop("If maxLPD is bigger than 0, it should be a whole number. Either define a number between 0 and 1 to use a proportion of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
    }
    maxLPD <- as.integer(maxLPD)
    if (maxLPD > n_samples) {
      stop("maxLPD cannot be bigger than the number of training samples. Either define a number between 0 and 1 to use a proportion of the number of training samples for the LPD calculation or a whole number larger than 1 and smaller than the number of training samples.")
    }
  }
  return(maxLPD)
}

################################################################################
# Categorial variable handling
################################################################################
.get_categorical_variables <- function(df, variables) {
  vars <- names(df)[which(sapply(df[, variables], class) %in% c("factor","character"))]
  return(vars)
}

.drop_unknown_levels <- function(reference, query = NULL, var) {
  reference[[var]] <- factor(droplevels(reference[[var]]))
  if (!is.null(query)) {
    query[[var]] <- factor(query[[var]], levels = levels(reference[[var]])) # will set unknown levels to NA
  }
  return(list(reference=reference, query=query))
}

.create_dummy_variables <- function(reference, query = NULL, weight = NULL, var) {
  # drop unknown levels and set unused levels to NA in reference and query
  result <- .drop_unknown_levels(reference, query, var)
  reference <- result$reference
  query <- result$query

  dvi_reference <- predict(caret::dummyVars(paste0("~", var), data = reference), reference)

  if (!is.null(query)) {
    dvi_query <- predict(caret::dummyVars(paste0("~", var), data = reference), query)
    dvi_query[is.na(query[ ,var]), ] <- 0
    query <- data.frame(query, dvi_query)
    query <- query[ , -which(names(query) == var)]
  }

  if (!is.null(weight)) {
    addweights <- data.frame(t(rep(weight[, which(names(weight)==var)], ncol(dvi_reference))))
    names(addweights) <- colnames(dvi_reference)
    weight <- data.frame(weight, addweights)
    weight <- weight[, -which(names(weight) == var)]
  }

  reference <- data.frame(reference, dvi_reference)
  reference <- reference[ , -which(names(reference) == var)]
  return(list(reference = reference, query = query, weight = weight))
}

.convert_factors_to_dummy <- function(reference, query = NULL, weight = NULL, variables = NULL) {
  if (is.null(variables) || length(variables) == 0) {
    return(list(reference = reference, query = query, weight = weight))
  }
  res <- list(reference = reference, query = query, weight = weight)
  for (var in variables) {
    res <- .create_dummy_variables(reference = res$reference, query = res$query, weight = res$weight, var = var)
  }
  return(res)
}

.prepare_categorical_variables <- function(reference, query = NULL, weight= NULL, variables) {
  vars <- .get_categorical_variables(reference, variables)
  res <- .convert_factors_to_dummy(reference, query, weight, vars)
  # we need to add new weights for the dummy variables
  return(list(reference = res$reference, query = res$query, weight = res$weight, catvars = vars))
}

################################################################################
# Weights handling
################################################################################
.prepare_weights <- function(weight = NULL, model = NULL, variables, useWeight = TRUE) {
  # default weights
  default_w <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
  names(default_w) <- variables

  if (!useWeight) {
    message("variable are not weighted. see ?aoa")
    return(default_w)
  }

  # fall back to defaults if weights are not provided
  if (is.null(weight)) {
     weight <- default_w
   }

  # model-based weights take precedence
  if (!is.null(model)) {
    weight <- .caret_get_weights(model, variables = variables)
  }

  weight <- .check_weights(weight, variables)
  return(weight)
}

# check if weights are correctly specified and set to 1 if not correctly specified
.check_weights <- function(weight, variables){
  if(inherits(weight, "list")){
    weight = as.data.frame(weight)
  }
  if (!all(variables %in% names(weight))) {
    stop("weights are not correctly specified. Check if all variables are included in the weight data.frame and if the column names match the variable names. See ?aoa")
  }
  if(nrow(weight) != 1) {
    message("variable weights are not correctly specified and will be ignored. See ?aoa")
    weight <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
    names(weight) <- variables
  }
  # subset weight to variables
  weight <- weight[, variables]

  if (any(weight < 0 )){
    weight[weight < 0] <- 0
    message("negative weights were set to 0")
  }
  if(sum(weight) == 0){
    stop("weights sum to 0. Check variable importance of the model, define weights manually or set useWeight=FALSE")
  }
  return(weight)
}

.apply_weights <- function(df, weight = NULL){
  if (is.null(weight)) {
    return(df)
  }
  if (nrow(weight) != 1) {
    stop("Expected weight to be a single row data frame.")
  }
  if (all(names(weight) %in% names(df))) {
    weight <- weight[, names(df), drop = FALSE]
  } else {
    stop("Weight columns do not match data frame columns.")
  }
  # multiply df cols by weight cols
  df <- sweep(df, 2, as.numeric(weight), FUN="*")
  return(df)
}


################################################################################
# Fold handling
################################################################################
.prepare_folds <- function(model = NULL, CVtrain = NULL, CVtest = NULL, useCV = FALSE) {
  if (!useCV) {
    message(
      "note: useCV is set to FALSE. The DI threshold is therefore based on all training data"
    )
    return(NULL)
  }

  if (!is.null(model)) {
    folds <- .caret_get_folds(model)
    return(folds)
  }

  if (is.null(CVtrain) && is.null(CVtest)) {
    message(
      "note: No model and no CV folds were given. The DI threshold is therefore based on all training data"
    )
    return(NULL)
  }

  has_CVtrain <- !is.null(CVtrain)
  has_CVtest <- !is.null(CVtest)

  # CVtest to list if it is a vector
  if (has_CVtest && !is.list(CVtest)) {
    CVtest <- split(seq_along(CVtest), CVtest)
  }

  if (!has_CVtest) {
    # difference of all samples and CVtrain
    CVtest <- lapply(CVtrain, function(x) {
      setdiff(seq_along(unlist(CVtrain)), x)
    })
  }
  if (!has_CVtrain) {
    # difference of all samples and CVtest
    CVtrain <- lapply(CVtest, function(x) setdiff(seq_along(unlist(CVtest)), x))
  }
  if (length(CVtrain) != length(CVtest)) {
    stop("CVtrain and CVtest should have the same number of folds")
  }
  folds <- list(CVtrain = CVtrain, CVtest = CVtest)
  return(folds)
}

################################################################################
# Variable handling
################################################################################
.prepare_variables <- function(variables, model, train) {

  if (length(variables) == 1) {  
    if (variables == "all") variables <- names(train)
    # takes precedence over train if model is given
    if (!is.null(model)) variables <- .caret_get_variables(model) 
  }

  if (!all(variables %in% names(train))) {
    stop("some of the specified variables are not included in the training data. Check variable names and if variables are included in the model. See ?aoa")
  }

  return(variables)
}


################################################################################
# Dissimlarity index calculation
################################################################################

# helper to derive DI threshold based on the distribution of DI in the training data
.di_threshold <- function(trainDI){
  threshold_quantile <- stats::quantile(trainDI, 0.75, na.rm = TRUE)
  threshold_iqr <- (1.5 * stats::IQR(trainDI, na.rm = TRUE))
  thres <- threshold_quantile + threshold_iqr
  # make sure that the threshold is not larger than the maximum DI in the training data
  thres <- min(thres, max(trainDI, na.rm = TRUE)) 
  return(thres)
}

