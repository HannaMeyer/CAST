#' Plot TrainDI
#' @param x trainDI object
#' @param ... other params
#'
#' @author Marvin Ludwig, Hanna Meyer
#' @export

plot.trainDI = function(x, ...){
  ggplot(data.frame(TrainDI = x$trainDI), aes_string(x = "TrainDI"))+
    geom_density()+
    geom_vline(xintercept = x$lower_threshold, linetype = "dashed")+
    geom_vline(xintercept = x$threshold, linetype = "dashed")+
    theme_bw()

}





#' Plot AOA
#'
#' @description Density plot of the DI with AOA threshold
#'
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
    targetDI = data.frame(DI = targetDI,
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

