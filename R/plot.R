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



#' @name plot
#' @param x An object of type \emph{nndm}.
#' @param ... other arguments.
#' @author Carles Milà
#'
#' @export
plot.nndm <- function(x, ...){

  # Prepare data for plotting: Gij function
  Gij_df <- data.frame(r=x$Gij[order(x$Gij)])
  Gij_df$val <- 1:nrow(Gij_df)/nrow(Gij_df)
  Gij_df <- Gij_df[Gij_df$r <= x$phi,]
  Gij_df <- rbind(Gij_df, data.frame(r=0, val=0))
  Gij_df <- rbind(Gij_df, data.frame(r=x$phi,
                                     val=sum(x$Gij<=x$phi)/length(x$Gij)))
  Gij_df$Function <- "1_Gij(r)"

  # Prepare data for plotting: Gjstar function
  Gjstar_df <- data.frame(r=x$Gjstar[order(x$Gjstar)])
  Gjstar_df$val <- 1:nrow(Gjstar_df)/nrow(Gjstar_df)
  Gjstar_df <- Gjstar_df[Gjstar_df$r <= x$phi,]
  Gjstar_df <- rbind(Gjstar_df, data.frame(r=0, val=0))
  Gjstar_df <- rbind(Gjstar_df, data.frame(r=x$phi,
                                           val=sum(x$Gjstar<=x$phi)/length(x$Gjstar)))
  Gjstar_df$Function <- "2_Gjstar(r)"

  # Prepare data for plotting: G function
  Gj_df <- data.frame(r=x$Gj[order(x$Gj)])
  Gj_df$val <- 1:nrow(Gj_df)/nrow(Gj_df)
  Gj_df <- Gj_df[Gj_df$r <= x$phi,]
  Gj_df <- rbind(Gj_df, data.frame(r=0, val=0))
  Gj_df <- rbind(Gj_df, data.frame(r=x$phi,
                                   val=sum(x$Gj<=x$phi)/length(x$Gj)))
  Gj_df$Function <- "3_Gj(r)"

  # Merge data for plotting
  Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)

  # Plot
  ggplot2::ggplot(Gplot) +
    ggplot2::geom_step(ggplot2::aes_string(x="r", y="val", colour="Function", size="Function"),
                       alpha = 0.8) +
    ggplot2::scale_size_manual(values=c(1.1, 1.1, 0.5),
                               labels=c(expression(hat(G)[ij](r)),
                                        expression(hat(G)[j]^"*"*"(r,"*bold(L)*")"),
                                        expression(hat(G)[j](r)))) +
    ggplot2::scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9"),
                                 labels=c(expression(hat(G)[ij](r)),
                                          expression(hat(G)[j]^"*"*"(r,"*bold(L)*")"),
                                          expression(hat(G)[j](r)))) +
    ggplot2::ylab(expression(paste(hat(G)[ij](r), ", ",
                                   hat(G)[j]^"*"*"(r,"*bold(L)*")", ", ",
                                   hat(G)[j](r)))) +
    ggplot2::labs(colour="", size="") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.text.align=0,
                   legend.text=ggplot2::element_text(size=12))
}

