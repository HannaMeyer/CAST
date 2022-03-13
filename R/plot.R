#' Plot CAST classes
#' @description Generic plot function for trainDI and aoa
#'
#' @name plot
#' @param x trainDI object
#' @param ... other params
#'
#'
#' @author Marvin Ludwig, Hanna Meyer
#' @export

plot.trainDI = function(x, ...){
  ggplot(data.frame(TrainDI = x$trainDI), aes_string(x = "TrainDI"))+
    geom_density()+
    geom_vline(aes(xintercept = x$threshold, linetype = "AOA_threshold"))+
    scale_linetype_manual(name = "", values = c(AOA_threshold = "dashed"))+
    theme_bw()+
    theme(legend.position="bottom")

}








#' @name plot
#'
#' @param x aoa object
#' @param samplesize numeric. How many prediction samples should be plotted?
#' @param ... other params
#'
#' @import ggplot2
#'
#' @author Marvin Ludwig, Hanna Meyer
#'
#' @export

plot.aoa = function(x, samplesize = 1000, ...){


  trainDI = data.frame(DI = x$parameters$trainDI,
                       what = "trainDI")



  if(inherits(x$AOA, "RasterLayer")){
    targetDI = raster::sampleRegular(x$DI, size = samplesize)
    targetDI = data.frame(DI = as.numeric(targetDI),
                          what = "predictionDI")
  }else if(inherits(x$AOA, "stars")){
    targetDI = raster::sampleRegular(methods::as(x$DI, "Raster"), size = samplesize)
    targetDI = data.frame(DI = as.numeric(targetDI),
                          what = "predictionDI")
  }else if(inherits(x$AOA, "SpatRaster")){
    targetDI = raster::sampleRegular(methods::as(x$DI, "Raster"), size = samplesize)
    targetDI = data.frame(DI = as.numeric(targetDI),
                          what = "predictionDI")
  }else{
    targetDI = data.frame(DI = sample(x$DI, size = samplesize),
                          what = "predictionDI")
  }



  dfDI = rbind(trainDI, targetDI)


  ggplot(dfDI, aes_string(x = "DI", group = "what", fill = "what"))+
    geom_density(adjust=1.5, alpha=.4)+
    scale_fill_discrete(name = "Set")+
    geom_vline(aes(xintercept = x$parameters$threshold, linetype = "AOA_threshold"))+
    scale_linetype_manual(name = "", values = c(AOA_threshold = "dashed"))+
    theme_bw()+
    theme(legend.position = "bottom")
}

