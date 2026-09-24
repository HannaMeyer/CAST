# Get trainingdata from train object
.caret_get_data <- function(model) {
  df <- as.data.frame(model$trainingData)
  return(df)
}

# Get variables from train object
.caret_get_variables <- function(model) {
  vars <- names(model$trainingData)
  if (is.null(vars)) {
    stop("no variable names could be retrieved from the given model. Check model$trainingData")
  }
  drop_idx <- which(vars == ".outcome")
  return(vars[-drop_idx])
}

.caret_get_folds <- function(model) {
  if (tolower(model$control$method) != "cv") {
    message(
      "note: No CV was used for model training. The DI threshold is therefore based on all training data"
    )
    return(list(CVtrain = NULL, CVtest = NULL))
  } else {
    CVtest <- model$control$indexOut
    CVtrain <- model$control$index
    return(list(CVtrain = CVtrain, CVtest = CVtest))
  }
}

# Get weights from train object
.caret_get_weights = function(model, variables){
  is_classification <- model$modelType == "Classification"
  if(is_classification){
    weight <- try(as.data.frame(t(apply(caret::varImp(model, scale=FALSE)$importance, 1, mean))))
  }else{
    weight <- try(as.data.frame(t(caret::varImp(model, scale=FALSE)$importance[, "Overall"])))
  }
  if (inherits(weight, "error") || ncol(weight) == 0) {
    weight <- as.data.frame(matrix(1, nrow = 1, ncol = length(variables)))
    names(weight) <- variables
    message("note: no variable importance could be retrieved from the given model. Variables are not weighted. Check caret::varImp(model)")
  } else {
    names(weight) <- rownames(caret::varImp(model, scale=FALSE)$importance)
  }
  return(.check_weights(weight, variables))
}


