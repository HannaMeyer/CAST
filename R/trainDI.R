#' Dissimilarity Index
#' @description
#' This function estimates the Dissimilarity Index (DI) of
#' of spatial prediction models. Predictors can be weighted based on the internal
#' variable importance of the machine learning algorithm used for model training.
#'
#' @param model A train object created with caret used to extract weights from (based on variable importance) as well as cross-validation folds
#' @param train A data.frame containing the data used for model training. Only required when no model is given
#' @param weight A data.frame containing weights for each variable. Only required if no model is given.
#' @param variables character vector of predictor variables. if "all" then all variables
#' of the model are used or if no model is given then of the train dataset.
#' @param folds Numeric or character. Folds for cross validation. E.g. Spatial cluster affiliation for each data point.
#' Should be used if replicates are present. Only required if no model is given.
#'
#' @seealso \code{\link{aoa}}
#' @importFrom graphics boxplot
#' @import ggplot2
#'
#' @return A trainDI object
#' @export trainDI
#'
#' @author
#' Hanna Meyer, Marvin Ludwig
#'
#' @references Meyer, H., Pebesma, E. (2020): Predicting into unknown space?
#' Estimating the area of applicability of spatial prediction models.
#' \url{https://arxiv.org/abs/2005.07939}


trainDI = function(model,
                   train = NULL,
                   variables = "all",
                   weight = NA,
                   folds = NULL){

  # get parameters if they are not provided in function call-----
  if(is.null(train)){train = aoa_get_train(model)}
  if(variables == "all"){variables = aoa_get_variables(variables, model, train)}
  if(is.na(weight)){weight = aoa_get_weights(model, variables = variables)}
  if(is.null(folds)){folds = aoa_get_folds(model, folds)}

  # check for input errors -----
  if(nrow(train)<=1){stop("at least two training points need to be specified")}

  # reduce train to specified variables
  train <- train[,na.omit(match(variables, names(train)))]

  train_backup = train

  # convert categorial variables
  catupdate = aoa_categorial_train(train, variables, weight)

  train = catupdate$train
  weight = catupdate$weight

  # scale train
  train <- scale(train)

  # make sure all variables have variance
  if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
    stop("some variables in train seem to have no variance")
  }

  # save scale param for later
  scaleparam <- attributes(train)


  # multiply train data with variable weights (from variable importance)
  if(!inherits(weight, "error")){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }


  # calculate average mean distance between training data

  ### Great potential for parallel computing -----

  trainDist_avrg = c()
  trainDist_min = c()

  for(i in seq(nrow(train))){

    # these two distance matrices should give the same output
    # it should be enough to calculate the second one, mask the self-distance as NA
    # (also important for fold masking afterward that the rows dont get mixed up)
    # and calculate the average afterwards


    # distance to all other training data (for average)
    trainDistAll = FNN::knnx.dist(train, t(train[i,]), k = nrow(train))[-1]
    trainDist_avrg = append(trainDist_avrg, mean(trainDistAll, na.rm = TRUE))

    # calculate distance to other training data:
    trainDist <- FNN::knnx.dist(t(matrix(train[i,])),train,k=1)
    trainDist[i] <- NA


    # mask of other folds
    if (!is.null(folds)){
      trainDist[folds==folds[i]] <- NA
    }

    trainDist_min = append(trainDist_min, min(trainDist, na.rm = TRUE))

  }
  trainDist_avrgmean = mean(trainDist_avrg,na.rm=TRUE)




  # Train Dissimilarity Index -----
  TrainDI <- trainDist_min/trainDist_avrgmean
  AOA_train_stats <- quantile(TrainDI,
                              probs = c(0.25,0.5,0.75,0.9,0.95,0.99,1),na.rm = TRUE)

  # AOA Threshold ----
  thres <- grDevices::boxplot.stats(TrainDI)$stats[5]
  lower_thres <- grDevices::boxplot.stats(TrainDI)$stats[1]

  # Return: AOA Object -------

  aoa_results = list(
    train = train_backup,
    weight = weight,
    variables = variables,
    catvars = catupdate$catvars,
    scaleparam = scaleparam,
    #trainDist_avrg = trainDist_avrg,
    trainDist_avrgmean = trainDist_avrgmean,
    #trainDist_min = trainDist_min,
    AOA_train_stats = AOA_train_stats,
    thres = thres,
    lower_thres = lower_thres,
    TrainDI = TrainDI
  )

  class(aoa_results) = "trainDI"

  return(aoa_results)
}

