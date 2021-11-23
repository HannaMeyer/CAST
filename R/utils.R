#' Encode categorial variables
#' @noRd


aoa_categorial_train <- function(train, variables, weight){

  # get all categorial variables
  catvars <- tryCatch(names(train)[which(sapply(train[,variables], class)%in%c("factor","character"))],
                      error=function(e) e)

  if (!inherits(catvars,"error")&length(catvars)>0){
    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      train[,catvar]<-droplevels(train[,catvar])

      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = train),train)
      train <- data.frame(train,dvi_train)

      if(!inherits(weight, "error")){
        addweights <- data.frame(t(rep(weight[,which(names(weight)==catvar)],
                                       ncol(dvi_train))))
        names(addweights)<- colnames(dvi_train)
        weight <- data.frame(weight,addweights)
      }
    }
    if(!inherits(weight, "error")){
      weight <- weight[,-which(names(weight)%in%catvars)]
    }
    train <- train[,-which(names(train)%in%catvars)]
  }
  return(list(train = train, weight = weight, catvars = catvars))


}



#' Get weights from train object
#' @noRd

aoa_get_weights = function(model, variables){

  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error")){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    # set all weights to 1
    weight <- as.data.frame(t(rep(1, length(variables))))
    names(weight) = variables
    message("note: variables were not weighted either because no weights or model were given,
    or no variable importance could be retrieved from the given model.
    Check caret::varImp(model)")
  }

  #set negative weights to 0
  if(!inherits(weight, "error")){
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  return(weight)

}


#' Get trainingdata from train object
#' @noRd

aoa_get_train <- function(model){

  train <- as.data.frame(model$trainingData)
  return(train)


}


#' Get folds from train object
#' @noRd

aoa_get_folds <- function(model, folds){
  if(!is.null(model)&is.null(folds)){
    CVfolds <- tryCatch(reshape::melt(model$control$indexOut),
                        error=function(e) e)
    if(inherits(CVfolds, "error")){
      message("note: Either no model was given or no CV was used for model training. The DI threshold is therefore based on all training data")
    }else{
      CVfolds <- CVfolds[order(CVfolds$value),]
    }
  }

  return(CVfolds$L1)

}


#' Get variables from train object
#' @noRd

aoa_get_variables <- function(variables, model, train){

  if(length(variables) == 1){
    if(variables == "all"){
      if(!is.na(model)[1]){
        variables <- names(model$trainingData)[-which(names(model$trainingData)==".outcome")]
      }else{
        variables <- names(train)
      }
    }
  }
  return(variables)


}
