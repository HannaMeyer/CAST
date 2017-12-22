#' Create Space-time Folds
#' @description Create spatial, temporal or spatio-temporal Folds for cross validation
#' @param x data.frame containing spatio-temporal data
#' @param spacevar Character indicating which column of x identifies the
#' spatial units (e.g. ID of weather stations)
#' @param timevar Character indicating which column of x identifies the
#' temporal units (e.g. the day of the year)
#' @param k numeric. Number of folds. If spacevar or timevar is NA and a
#' leave one location out or leave one time step out cv should be performed,
#' set k to the number of unique spatial or temporal units.
#' @param seed numeric. See ?seed
#' @return A list that contains a list for model training and a list for
#' model validation that can directly be used as "index" and "indexOut" in
#' caret's train function
#' @author Hanna Meyer
#' @examples
#' library(GSIF)
#' data(cookfarm)
#' indices <- CreateSpacetimeFolds(cookfarm$readings,"SOURCEID","Date")
#' str(indices)
#' ### Prepare for leave one location out cross validation
#' indices <- CreateSpacetimeFolds(cookfarm$readings,spacevar="SOURCEID",
#' k=length(unique(cookfarm$readings$SOURCEID)))
#' @export CreateSpacetimeFolds
#' @aliases CreateSpacetimeFolds

CreateSpacetimeFolds <- function(x,spacevar=NA,timevar=NA,
                                 k=10,seed=sample(1:1000, 1)){
  #split space into k folds
  if(!is.na(spacevar)){
    set.seed(seed)
  spacefolds <- lapply(caret::createFolds(1:length(unique(x[,spacevar])),k),function(y){
    unique(x[,spacevar])[y]})
  }
  #split time into k folds
  if(!is.na(timevar)){
    set.seed(seed)
  timefolds <- lapply(caret::createFolds(1:length(unique(x[,timevar])),k),function(y){
    unique(x[,timevar])[y]})
  }
  # combine space and time folds
  cvindices_train <- list()
  cvindices_test <- list()
  for (i in 1:k){
    if(!is.na(timevar)&!is.na(spacevar)){
    cvindices_test[[i]]<- which(x[,spacevar]%in%spacefolds[[i]]&
                                  x[,timevar]%in%timefolds[[i]])
    cvindices_train[[i]]<- which(!x[,spacevar]%in%spacefolds[[i]]&
                                   !x[,timevar]%in%timefolds[[i]])
    }
    if(is.na(timevar)&!is.na(spacevar)){
      cvindices_test[[i]]<- which(x[,spacevar]%in%spacefolds[[i]])
      cvindices_train[[i]]<- which(!x[,spacevar]%in%spacefolds[[i]])
    }
    if(!is.na(timevar)&is.na(spacevar)){
      cvindices_test[[i]]<- which(x[,timevar]%in%timefolds[[i]])
      cvindices_train[[i]]<- which(!x[,timevar]%in%timefolds[[i]])
    }
  }
  return(list("index"=cvindices_train,"indexOut"=cvindices_test))
}
