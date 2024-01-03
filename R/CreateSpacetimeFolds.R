#' Create Space-time Folds
#' @description Create spatial, temporal or spatio-temporal Folds for cross validation based on pre-defined groups
#' @param x data.frame containing spatio-temporal data
#' @param spacevar Character indicating which column of x identifies the
#' spatial units (e.g. ID of weather stations)
#' @param timevar Character indicating which column of x identifies the
#' temporal units (e.g. the day of the year)
#' @param k numeric. Number of folds. If spacevar or timevar is NA and a
#' leave one location out or leave one time step out cv should be performed,
#' set k to the number of unique spatial or temporal units.
#' @param class Character indicating which column of x identifies a class unit (e.g. land cover)
#' @param seed numeric. See ?seed
#' @return A list that contains a list for model training and a list for
#' model validation that can directly be used as "index" and "indexOut" in
#' caret's trainControl function
#' @details The function creates train and test sets by taking (spatial and/or temporal) groups into account.
#' In contrast to \code{\link{nndm}}, it requires that the groups are already defined (e.g. spatial clusters or blocks or temporal units).
#' Using "class" is helpful in the case that data are clustered in space
#' and are categorical. E.g This is the case for land cover classifications when
#' training data come as training polygons. In this case the data should be split in a way
#' that entire polygons are held back (spacevar="polygonID") but at the same time the distribution of classes
#' should be similar in each fold (class="LUC").
#' @note Standard k-fold cross-validation can lead to considerable misinterpretation in spatial-temporal modelling tasks. This function can be used to prepare a Leave-Location-Out, Leave-Time-Out or Leave-Location-and-Time-Out cross-validation as target-oriented validation strategies for spatial-temporal prediction tasks. See Meyer et al. (2018) for further information.
#' @author Hanna Meyer
#' @seealso \code{\link{trainControl}},\code{\link{ffs}}, \code{\link{nndm}}
#' @references
#' Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018): Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software 101: 1-9.
#' @examples
#' \dontrun{
#' dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
#' ### Prepare for 10-fold Leave-Location-and-Time-Out cross validation
#' indices <- CreateSpacetimeFolds(dat,"SOURCEID","Date")
#' str(indices)
#' ### Prepare for 10-fold Leave-Location-Out cross validation
#' indices <- CreateSpacetimeFolds(dat,spacevar="SOURCEID")
#' str(indices)
#' ### Prepare for leave-One-Location-Out cross validation
#' indices <- CreateSpacetimeFolds(dat,spacevar="SOURCEID",
#'     k=length(unique(dat$SOURCEID)))
#' str(indices)
#' }
#' @export CreateSpacetimeFolds
#' @aliases CreateSpacetimeFolds

CreateSpacetimeFolds <- function(x,spacevar=NA,timevar=NA,
                                 k=10,class=NA,seed=sample(1:1000, 1)){
  x <- data.frame(x)
  ### if classification is used, make sure that classes are equally distributed across folds
  if(!is.na(class)){
    unit <- unique(x[,c(spacevar,class)])
    unit$CAST_fold <- createFolds(unit[,which(names(unit)==class)],k = k,list=FALSE)
    #x <- merge(x,unit,by.x=c(spacevar,class),by.y=c(spacevar,class),all.x=TRUE,sort=FALSE)
    x <- plyr::join(x,unit,by=c(spacevar,class),match="first")
    spacevar <- "CAST_fold"  }

  if(!is.na(spacevar)){
    if(k>length(unique(x[,spacevar]))){
      k <- length(unique(x[,spacevar]))
      print(paste0("warning: k is higher than number of unique locations. k is set to ",k))
    }
  }
  if(!is.na(timevar)){
    if(k>length(unique(x[,timevar]))){
      k <- length(unique(x[,timevar]))
      print(paste0("warning: k is higher than number of unique points in time. k is set to ",k))
    }
  }
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
