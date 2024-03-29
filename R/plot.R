#' Plot CAST classes
#' @description Generic plot function for CAST Classes
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
#' @param variable character. Variable for which to generate the density plot. 'DI' or 'LPD'
#' @param ... other params
#'
#' @import ggplot2
#'
#' @author Marvin Ludwig, Hanna Meyer
#'
#' @export

plot.aoa = function(x, samplesize = 1000, variable = "DI", ...){

  if (variable == "DI") {
    trainDI = data.frame(DI = x$parameters$trainDI,
                         what = "trainDI")

    if(inherits(x$AOA, "RasterLayer")){
      targetDI = terra::spatSample(methods::as(x$DI, "SpatRaster"),
                                   size = samplesize, method = "regular")
      targetDI = data.frame(DI = as.numeric(targetDI[, 1]),
                            what = "predictionDI")
    }else if(inherits(x$AOA, "stars")){
      targetDI = terra::spatSample(methods::as(x$DI, "SpatRaster"),
                                   size = samplesize, method = "regular")
      targetDI = data.frame(DI = as.numeric(targetDI[, 1]),
                            what = "predictionDI")
    }else if(inherits(x$AOA, "SpatRaster")){
      targetDI = terra::spatSample(x$DI, size = samplesize, method = "regular")
      targetDI = data.frame(DI = as.numeric(targetDI[, 1]),
                            what = "predictionDI")
    }else{
      targetDI = data.frame(DI = sample(x$DI, size = samplesize),
                            what = "predictionDI")
    }

    dfDI = rbind(trainDI, targetDI)

    plot = ggplot(dfDI, aes_string(x = "DI", group = "what", fill = "what"))+
      geom_density(adjust=1.5, alpha=.4)+
      scale_fill_discrete(name = "Set")+
      geom_vline(aes(xintercept = x$parameters$threshold, linetype = "AOA_threshold"))+
      scale_linetype_manual(name = "", values = c(AOA_threshold = "dashed"))+
      theme_bw()+
      theme(legend.position = "bottom")
  }


  if (variable == "LPD") {
    trainLPD = data.frame(LPD = x$parameters$trainLPD,
                          what = "trainLPD")



    if(inherits(x$AOA, "RasterLayer")){
      targetLPD = terra::spatSample(methods::as(x$LPD, "SpatRaster"),
                                    size = samplesize, method = "regular")
      targetLPD = data.frame(LPD = as.numeric(targetLPD[, 1]),
                             what = "predictionLPD")
    }else if(inherits(x$AOA, "stars")){
      targetLPD = terra::spatSample(methods::as(x$LPD, "SpatRaster"),
                                    size = samplesize, method = "regular")
      targetLPD = data.frame(LPD = as.numeric(targetLPD[, 1]),
                             what = "predictionLPD")
    }else if(inherits(x$AOA, "SpatRaster")){
      targetLPD = terra::spatSample(x$LPD, size = samplesize, method = "regular")
      targetLPD = data.frame(LPD = as.numeric(targetLPD[, 1]),
                             what = "predictionLPD")
    }else{
      targetLPD = data.frame(LPD = sample(x$LPD, size = samplesize),
                             what = "predictionLPD")
    }

    dfLPD = rbind(trainLPD, targetLPD)


    plot = ggplot(dfLPD, aes_string(x = "LPD", group = "what", fill = "what"))+
      geom_density(adjust=1.5, alpha=0.4)+
      scale_fill_discrete(name = "Set")+
      geom_vline(aes(xintercept = median(x$parameters$trainLPD), linetype = "MtrainLPD"))+
      scale_linetype_manual(name = "", values = c(MtrainLPD = "dashed"))+
      theme_bw()+
      theme(legend.position = "bottom")

  }

  return(plot)
}



#' @name plot
#' @param x An object of type \emph{nndm}.
#' @param type String, defaults to "strict" to show the original nearest neighbour distance definitions in the legend.
#' Alternatively, set to "simple" to have more intuitive labels.
#' @param ... other arguments.
#' @author Carles Milà
#'
#' @export
plot.nndm <- function(x, type="strict", stat = "ecdf", ...){

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

  # Merge data for plotting, get maxdist relevant for plotting
  if(any(Gj_df$val==1)&any(Gjstar_df$val==1)&any(Gij_df$val==1)){
    Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)
    maxdist <- max(Gplot$r[Gplot$val!=1]) + 1e-9
    Gplot <- Gplot[Gplot$r <= maxdist,]
    Gplot <- rbind(Gplot, data.frame(r=maxdist, val=1,
                                     Function = c("1_Gij(r)", "2_Gjstar(r)", "3_Gj(r)")))
  }else{
    Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)
  }

  # Define colours matching those of geodist
  myColors <- RColorBrewer::brewer.pal(3, "Dark2")

  # Plot
  if(stat=="ecdf"){
    p <- ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", col="Function")) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step", lwd = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::ylab("ECDF") +
      ggplot2::labs(group="Distance function", col="Distance function") +
      ggplot2::theme(legend.position = "bottom",
                     legend.text=ggplot2::element_text(size=10))

    if(type=="strict"){
      p <-  p +
        ggplot2::scale_colour_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                     labels=c(expression(hat(G)[ij](r)),
                                              expression(hat(G)[j]^"*"*"(r,L)"),
                                              expression(hat(G)[j](r))))
    }else if(type == "simple"){
      p <-  p +
        ggplot2::scale_colour_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                     labels=c("prediction-to-sample",
                                              "CV-distances",
                                              "sample-to-sample"))
    }

  }else if(stat=="density"){
    p <- ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", fill="Function")) +
      ggplot2::geom_density(adjust=1.5, alpha=.5, stat=stat, lwd = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Density") +
      ggplot2::labs(group="Distance function", col="Distance function") +
      ggplot2::theme(legend.position = "bottom",
                     legend.text=ggplot2::element_text(size=10))

    if(type=="strict"){
      p <-  p +
        ggplot2::scale_fill_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                   labels=c(expression(hat(G)[ij](r)),
                                            expression(hat(G)[j]^"*"*"(r,L)"),
                                            expression(hat(G)[j](r))))
    }else if(type == "simple"){
      p <-  p +
        ggplot2::scale_fill_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                   labels=c("prediction-to-sample",
                                            "CV-distances",
                                            "sample-to-sample"))
    }
  }

  p

}


#' @name plot
#' @param x An object of type \emph{knndm}.
#' @param type String, defaults to "strict" to show the original nearest neighbour distance definitions in the legend.
#' Alternatively, set to "simple" to have more intuitive labels.
#' @param stat String, defaults to "ecdf" but can be set to "density" to estimate density functions.
#' @param ... other arguments.
#' @author Carles Milà
#'
#' @export
plot.knndm <- function(x, type="strict", stat = "ecdf", ...){

  # Prepare data for plotting: Gij function
  Gij_df <- data.frame(r=x$Gij[order(x$Gij)])
  Gij_df$Function <- "1_Gij(r)"

  # Prepare data for plotting: Gjstar function
  Gjstar_df <- data.frame(r=x$Gjstar[order(x$Gjstar)])
  Gjstar_df$Function <- "2_Gjstar(r)"

  # Prepare data for plotting: G function
  Gj_df <- data.frame(r=x$Gj[order(x$Gj)])
  Gj_df$Function <- "3_Gj(r)"

  # Merge data for plotting
  Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)

  # Define colours matching those of geodist
  myColors <- RColorBrewer::brewer.pal(3, "Dark2")

  # Plot
  if(stat=="ecdf"){
    p <- ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", col="Function")) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step", lwd = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::ylab("ECDF") +
      ggplot2::labs(group="Distance function", col="Distance function") +
      ggplot2::theme(legend.position = "bottom",
                     legend.text=ggplot2::element_text(size=10))

    if(type=="strict"){
      p <-  p +
        ggplot2::scale_colour_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                     labels=c(expression(hat(G)[ij](r)),
                                              expression(hat(G)[j]^"*"*"(r,L)"),
                                              expression(hat(G)[j](r))))
    }else if(type == "simple"){
      p <-  p +
        ggplot2::scale_colour_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                     labels=c("prediction-to-sample",
                                              "CV-distances",
                                              "sample-to-sample"))
    }

  }else if(stat=="density"){
    p <- ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", fill="Function")) +
      ggplot2::geom_density(adjust=1.5, alpha=.5, stat=stat, lwd = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Density") +
      ggplot2::labs(group="Distance function", col="Distance function") +
      ggplot2::theme(legend.position = "bottom",
                     legend.text=ggplot2::element_text(size=10))

    if(type=="strict"){
      p <-  p +
        ggplot2::scale_fill_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                   labels=c(expression(hat(G)[ij](r)),
                                            expression(hat(G)[j]^"*"*"(r,L)"),
                                            expression(hat(G)[j](r))))
    }else if(type == "simple"){
      p <-  p +
        ggplot2::scale_fill_manual(values=c(myColors[2], myColors[3], myColors[1]),
                                   labels=c("prediction-to-sample",
                                            "CV-distances",
                                            "sample-to-sample"))
    }
  }

  p
}

#' Plot results of a Forward feature selection or best subset selection
#' @description A plotting function for a forward feature selection result.
#' Each point is the mean performance of a model run. Error bars represent
#' the standard errors from cross validation.
#' Marked points show the best model from each number of variables until a further variable
#' could not improve the results.
#' If type=="selected", the contribution of the selected variables to the model
#' performance is shown.
#' @param x Result of a forward feature selection see \code{\link{ffs}}
#' @param plotType character. Either "all" or "selected"
#' @param palette A color palette
#' @param reverse Character. Should the palette be reversed?
#' @param marker Character. Color to mark the best models
#' @param size Numeric. Size of the points
#' @param lwd Numeric. Width of the error bars
#' @param pch Numeric. Type of point marking the best models
#' @param ... Further arguments for base plot if type="selected"
#' @author Marvin Ludwig, Hanna Meyer
#' @examples
#' \dontrun{
#' data(splotdata)
#' splotdata <- st_drop_geometry(splotdata)
#' ffsmodel <- ffs(splotdata[,6:16], splotdata$Species_richness, ntree = 10)
#' plot(ffsmodel)
#' #plot performance of selected variables only:
#' plot(ffsmodel,plotType="selected")
#'}
#' @name plot
#' @importFrom forcats fct_rev fct_inorder
#' @export



plot.ffs <- function(x,plotType="all",palette=rainbow,reverse=FALSE,
                     marker="black",size=1.5,lwd=0.5,
                     pch=21,...){
  metric <- x$metric
  if (is.null(x$type)){
    x$type <- "ffs"
  }
  if(is.null(x$minVar)){
    x$minVar <- 2
  }
  if(x$type=="bss"&plotType=="selected"){
    type <- "all"
    print("warning: type must be 'all' for a bss model")
  }
  if (plotType=="selected"){

    plot_df = data.frame(labels = forcats::fct_rev(forcats::fct_inorder(c(paste(x$selectedvars[1:x$minVar], collapse = "\n + "),
                                    paste("+", x$selectedvars[-1:-x$minVar], sep = " ")))),
                         perf = x$selectedvars_perf,
                         perfse = x$selectedvars_perf_SE)


    p <- ggplot(plot_df, aes_string(x = "perf", y = "labels"))+
      geom_point()+
      geom_segment(aes_string(x = "perf - perfse", xend = "perf + perfse",
                       y = "labels", yend = "labels"))+
      xlab(x$metric)+
      ylab(NULL)
    return(p)

  }else{


    output_df <- x$perf_all
    output_df$run <- seq(nrow(output_df))
    names(output_df)[which(names(output_df)==metric)] <- "value"

    if (x$type=="bss"){
      bestmodels <- output_df$run[which(output_df$value==x$selectedvars_perf)]
    }else{

      bestmodels <- c()
      for (i in unique(output_df$nvar)){
        if (x$maximize){
          bestmodels <- c(bestmodels,
                          output_df$run[output_df$nvar==i][which(output_df$value[
                            output_df$nvar==i]==max(output_df$value[output_df$nvar==i]))][1])
        }else{
          bestmodels <- c(bestmodels,
                          output_df$run[output_df$nvar==i][which(output_df$value[
                            output_df$nvar==i]==min(output_df$value[output_df$nvar==i]))][1])
        }
      }
      bestmodels <- bestmodels[1:(length(x$selectedvars)-1)]
    }

    if (!reverse){
      cols <- palette(max(output_df$nvar)-(min(output_df$nvar)-1))
    }else{
      cols <- rev(palette(max(output_df$nvar)-(min(output_df$nvar)-1)))
    }
    ymin <- output_df$value - output_df$SE
    ymax <- output_df$value + output_df$SE
    if (max(output_df$nvar)>11){
      p <- ggplot2::ggplot(output_df, ggplot2::aes_string(x = "run", y = "value"))+
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                               color = cols[output_df$nvar-(min(output_df$nvar)-1)],lwd=lwd)+
        ggplot2::geom_point(ggplot2::aes_string(colour="nvar"),size=size)+
        ggplot2::geom_point(data=output_df[bestmodels, ],
                            ggplot2::aes_string(x = "run", y = "value"),
                            pch=pch,colour=marker,size=size)+
        ggplot2::scale_x_continuous(name = "Model run", breaks = pretty(output_df$run))+
        ggplot2::scale_y_continuous(name = metric)+
        ggplot2::scale_colour_gradientn(breaks=seq(2,max(output_df$nvar),
                                                   by=ceiling(max(output_df$nvar)/5)),
                                        colours = cols, name = "variables",guide = "colourbar")
    }else{
      dfint <- output_df
      dfint$nvar <- as.factor(dfint$nvar)
      p <- ggplot2::ggplot(dfint, ggplot2::aes_string(x = "run", y = "value"))+
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                               color = cols[output_df$nvar-(min(output_df$nvar)-1)],lwd=lwd)+
        ggplot2::geom_point(ggplot2::aes_string(colour="nvar"),size=size)+
        ggplot2::geom_point(data=output_df[bestmodels, ],
                            ggplot2::aes_string(x = "run", y = "value"),
                            pch=pch,colour=marker,size=size)+
        ggplot2::scale_x_continuous(name = "Model run", breaks = pretty(dfint$run))+
        ggplot2::scale_y_continuous(name = metric)+
        ggplot2::scale_colour_manual(values = cols, name = "variables")

    }

    return(p)
  }
}


#' @name plot
#' @description Density plot of nearest neighbor distances in geographic space or feature space between training data as well as between training data and
#' prediction locations.
#' Optional, the nearest neighbor distances between training data and test data or between training data and CV iterations is shown.
#' The plot can be used to check the suitability of a chosen CV method to be representative to estimate map accuracy.
#' @param x geodist, see \code{\link{geodist}}
#' @param unit character. Only if type=="geo" and only applied to the plot. Supported: "m" or "km".
#' @param stat "density" for density plot or "ecdf" for empirical cumulative distribution function plot.
#' @export
#' @return a ggplot
#'



plot.geodist <- function(x, unit = "m", stat = "density", ...){

  # Define colours - they must match those of knndm and nndm
  labs <- c("sample-to-sample",
            "prediction-to-sample",
            "CV-distances",
            "test-to-sample")
  myColors <- RColorBrewer::brewer.pal(4, "Dark2")
  names(myColors) <- labs


  type <- attr(x, "type")

  if(unit=="km"){
    x$dist <- x$dist/1000
    xlabs <- "geographic distances (km)"
  }else{
    xlabs <- "geographic distances (m)"
  }

  if( type=="feature"){ xlabs <- "feature space distances"}
  what <- "" #just to avoid check note
  if (type=="feature"){unit ="unitless"}

  if (type=="time"){unit = attr(x,"unit")}
  if( type=="time"){ xlabs <- paste0("temporal distances (",unit,")")}
  if(stat=="density"){
    p <- ggplot2::ggplot(data=x, aes(x=dist, group=what, fill=what)) +
      ggplot2::geom_density(adjust=1.5, alpha=.5, stat=stat, lwd = 0.3) +
      ggplot2::scale_fill_manual(name = "distance function", values = myColors) +
      ggplot2::theme_bw() +
      ggplot2::xlab(xlabs) +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }else if(stat=="ecdf"){
    p <- ggplot2::ggplot(data=x, aes(x=dist, group=what, col=what)) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step", lwd = 1) +
      ggplot2::scale_color_manual(name = "distance function", values = myColors) +
      ggplot2::theme_bw() +
      ggplot2::xlab(xlabs) +
      ggplot2::ylab("ECDF") +
      ggplot2::theme(legend.position="bottom",
                     plot.margin = unit(c(0,0.5,0,0),"cm"))
  }
  p
}


#' @name plot
#' @description Plot the DI/LPD and errormetric from Cross-Validation with the modeled relationship
#' @param x errorModel, see \code{\link{DItoErrormetric}}
#' @param ... other params
#' @export
#' @return a ggplot
#'


plot.errorModel <- function(x, ...){

  variable = attr(x, "variable")
  metric = attr(x, "metric")

  performance = attr(x, "performance")[,c(variable, "metric")]
  performance$what = "cross-validation"

  model_line = data.frame(variable = performance[, variable],
                          metric = predict(x, performance),
                          what = "model")

  p = ggplot()+
    geom_point(data = performance, mapping = aes_string(x = variable, y = "metric", shape = "what"))+
    geom_line(data = model_line, mapping =  aes_string(x = "variable", y = "metric", linetype = "what"), lwd = 1)+
    labs(x = variable, y = metric)+
    theme(legend.title = element_blank(), legend.position = "bottom")

  return(p)

}
