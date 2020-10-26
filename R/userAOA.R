#' User defined Area of Applicability
#' @description
#' @param AOA the result of ?aoa
#' @param model the model used to get the AOA
#' @param proficiency numeric. error accepted for the AOA. To first get an overview leave it to NA and select it after inspecting the plot
#' @param showplot Logical. Show visualization oder not?
#' @details
#' @return
#' @author
#' Hanna Meyer
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(caret)
#' library(viridis)
#' library(latticeExtra)
#'
#' # prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
#' pts <- st_as_sf(dat,coords=c("Easting","Northing"))
#' pts$ID <- 1:nrow(pts)
#' set.seed(100)
#' pts <- pts[1:30,]
#' studyArea <- stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))[[1:8]]
#' trainDat <- extract(studyArea,pts,df=TRUE)
#' trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
#'
#' # train a model:
#' set.seed(100)
#' variables <- c("DEM","NDRE.Sd","TWI")
#' model <- train(trainDat[,which(names(trainDat)%in%variables)],
#' trainDat$VW,method="rf",importance=TRUE,tuneLength=1,
#' trControl=trainControl(method="cv",number=5,savePredictions=TRUE))
#'
#' #...then calculate the AOA of the trained model for the study area:
#' AOA <- aoa(studyArea,model)
#' #...or get the area for a user-defined performance
#' AOA_new <- userAOA(AOA,model,proficiency=0.7)
#' }
#' @export userAOA
#' @aliases userAOA

userAOA <- function(AOA,model,proficiency=NULL,showplot=TRUE){
  for (i in 1:length(model$bestTune)){
    tunevar <- names(model$bestTune[i])
    preds <- model$pred[model$pred[,tunevar]==model$bestTune[,tunevar],] ###HIER FÜR ALLE TUNEVARS
    preds <- preds[order(preds$rowIndex),]
  }

  if(is.null(preds)){message("no cross-predictions can be retrieved from the model. Train with savePredictions=TRUE")}
  absError <- abs(preds$pred-preds$obs)
  th_orig <- attributes(AOA)$aoa_stats$threshold
  TrainDI <- attributes(AOA)$TrainDI$DI
  errormodel <- lm(absError~TrainDI)

  if(is.null(proficiency)){
    proficiency <-  predict(errormodel,data.frame("TrainDI"=th_orig))
  }
  AOAthres <- (-coefficients(summary(errormodel))[1]+proficiency)/coefficients(summary(errormodel))[2]

  AOA$AOA[AOA$DI>AOAthres] <- 0
  AOA$AOA[AOA$DI<=AOAthres] <- 1

  attributes(AOA)$aoa_stats$threshold <- AOAthres

  if(AOAthres>max(TrainDI)){message("Warning: new threshold > maximum DI of the training data. Assumptions about the threshold are made based on extrapolation.")}



    dat <- data.frame(TrainDI,absError)

    plt <- ggplot()+geom_point(aes(TrainDI,absError,colour="Training data"))+
      scale_color_manual(name="",values = c("Training data" = 'black'))+
      xlab("DI")+ylab("absolute error")+
      geom_abline(intercept = coefficients(summary(errormodel))[1],
                  slope =coefficients(summary(errormodel))[2],col="blue")+
      geom_vline(aes(xintercept=th_orig,linetype="original threshold"),col="red")+
      geom_vline(aes(xintercept=AOAthres,linetype="new threshold"),col="red")+
      xlim(c(min(TrainDI,AOAthres,na.rm=T), max(TrainDI,AOAthres,na.rm=T)))+
      ylim(c(min(absError,proficiency,na.rm=T),max(absError,proficiency,na.rm=T)))+
      scale_linetype_manual(name="",values=c("original threshold"="dashed",
                                             "new threshold"="solid"))
    if(showplot){
    print(plt)
    }
    return(list(AOA = AOA,
                plot = plt))
}

