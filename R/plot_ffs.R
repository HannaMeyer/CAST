#' Plot results of a Forward feature selection
#' @description A simple plotting function for a forward feature selection result
#' @param ffs_model Result of a forward feature selection see \code{\link{ffs}}
#' @param palette A color palette
#' @return A plot object. Grey bars represent standard errors from cross validation.
#' Marked points show the best model from each number of variables.
#' @author Marvin Ludwig and Hanna Meyer
#' @seealso \code{\link{ffs}}
#' @examples
#' \dontrun{
#' data(iris)
#' ffsmodel <- ffs(iris[,1:4],iris$Species)
#' plot_ffs(ffsmodel)
#'}
#' @export plot_ffs
#' @aliases plot_ffs



plot_ffs <- function(ffs_model,palette="viridis"){
  metric <- ffs_model$metric
  output_df <- ffs_model$perf_all
  output_df$run <- seq(nrow(output_df))
  names(output_df)[which(names(output_df)==metric)] <- "value"
  bestmodels <- c()
  for (i in unique(output_df$nvar)){
    if (ffs_model$maximize){
    bestmodels <- c(bestmodels,
                    output_df$run[output_df$nvar==i][which(output_df$value[
                      output_df$nvar==i]==max(output_df$value[output_df$nvar==i]))])
    }else{
      bestmodels <- c(bestmodels,
                      output_df$run[output_df$nvar==i][which(output_df$value[
                        output_df$nvar==i]==min(output_df$value[output_df$nvar==i]))])
    }
  }
  bestmodels <- bestmodels[1:(length(ffs_model$selectedvars)-1)]
  cols1 <- rev(eval(parse(text=paste0(palette,"(",max(output_df$nvar),")"))))
  cols <- cols1[output_df$nvar]
  output_df$ymin <- output_df$value - output_df$SE
  output_df$ymax <- output_df$value + output_df$SE
  col_border <- rep(NA,nrow(output_df))
  col_border[bestmodels] <- "black"


  ggplot2::ggplot(output_df, ggplot2::aes_string(x = "run", y = "value"))+
    ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "ymin", ymax = "ymax"), color = "gray90")+
    ggplot2::geom_point(ggplot2::aes_string(colour="nvar"))+
    ggplot2::geom_point(data=output_df[bestmodels, ],ggplot2::aes_string(x = "run", y = "value"),
                        pch=21,lwd=2)+
    ggplot2::scale_colour_gradientn(colours = cols1,name = "# of variables")+
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank())+
    ggplot2::scale_x_continuous(name = "Model run")+
    ggplot2::scale_y_continuous(name = metric)
}
