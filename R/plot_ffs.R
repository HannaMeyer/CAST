#' Plot results of a Forward feature selection
#' @description A plotting function for a forward feature selection result.
#' Each point is the mean performance of a model run. Error bars represent
#' the standard errors from cross validation.
#' Marked points show the best model from each number of variables until a further variable
#' could not improve the results.
#' @param ffs_model Result of a forward feature selection see \code{\link{ffs}}
#' @param palette A color palette
#' @param reverse Character. Should the palette be reversed?
#' @param marker Character. Color to mark the best models
#' @param size Numeric. Size of the points
#' @param lwd Numeric. Width of the error bars
#' @param pch Numeric. Type of point marking the best models
#' @author Marvin Ludwig and Hanna Meyer
#' @seealso \code{\link{ffs}}
#' @examples
#' \dontrun{
#' data(iris)
#' ffsmodel <- ffs(iris[,1:4],iris$Species)
#' plot_ffs(ffsmodel)
#'}
#' @export plot_ffs
#' @aliases plot_ffs plot_bss


plot_ffs <- function(ffs_model,palette=rainbow,reverse=FALSE, marker="black",size=1.5,lwd=0.5,
                     pch=21){
  metric <- ffs_model$metric
  output_df <- ffs_model$perf_all
  output_df$run <- seq(nrow(output_df))
  names(output_df)[which(names(output_df)==metric)] <- "value"

  if (ffs_model$type=="bss"){
    bestmodels <- output_df$run[which(output_df$value==ffs_model$selectedvars_perf)]
  }else{

  bestmodels <- c()
  for (i in unique(output_df$nvar)){
    if (ffs_model$maximize){
      bestmodels <- c(bestmodels,
                      output_df$run[output_df$nvar==i][which(output_df$value[
                        output_df$nvar==i]==max(output_df$value[output_df$nvar==i]))][1])
    }else{
      bestmodels <- c(bestmodels,
                      output_df$run[output_df$nvar==i][which(output_df$value[
                        output_df$nvar==i]==min(output_df$value[output_df$nvar==i]))][1])
    }
  }
  bestmodels <- bestmodels[1:(length(ffs_model$selectedvars)-1)]
}

  if (!reverse){
    cols <- palette(max(output_df$nvar)-1)
  }else{
    cols <- rev(palette(max(output_df$nvar)-1))
  }
  ymin <- output_df$value - output_df$SE
  ymax <- output_df$value + output_df$SE
  if (max(output_df$nvar)>11){
    p <- ggplot2::ggplot(output_df, ggplot2::aes_string(x = "run", y = "value"))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                             color = cols[output_df$nvar-1],lwd=lwd)+
      ggplot2::geom_point(ggplot2::aes_string(colour="nvar"),size=size)+
      ggplot2::geom_point(data=output_df[bestmodels, ],
                          ggplot2::aes_string(x = "run", y = "value"),
                          pch=pch,colour=marker,lwd=size)+
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
                             color = cols[output_df$nvar-1],lwd=lwd)+
      ggplot2::geom_point(ggplot2::aes_string(colour="nvar"),size=size)+
      ggplot2::geom_point(data=output_df[bestmodels, ],
                          ggplot2::aes_string(x = "run", y = "value"),
                          pch=pch,colour=marker,lwd=size)+
      ggplot2::scale_x_continuous(name = "Model run", breaks = pretty(dfint$run))+
      ggplot2::scale_y_continuous(name = metric)+
      ggplot2::scale_colour_manual(values = cols, name = "variables")

  }
  return(p)
}
