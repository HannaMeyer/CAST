#' Plot results of a Forward feature selection or best subset selection
#' @description plot_ffs() is deprecated and will be removed soon. Please use generic plot() function on ffs object.
#' A plotting function for a forward feature selection result.
#' Each point is the mean performance of a model run. Error bars represent
#' the standard errors from cross validation.
#' Marked points show the best model from each number of variables until a further variable
#' could not improve the results.
#' If type=="selected", the contribution of the selected variables to the model
#' performance is shown.
#' @param ffs_model Result of a forward feature selection see \code{\link{ffs}}
#' @param plotType character. Either "all" or "selected"
#' @param palette A color palette
#' @param reverse Character. Should the palette be reversed?
#' @param marker Character. Color to mark the best models
#' @param size Numeric. Size of the points
#' @param lwd Numeric. Width of the error bars
#' @param pch Numeric. Type of point marking the best models
#' @param ... Further arguments for base plot if type="selected"
#' @author Marvin Ludwig and Hanna Meyer
#' @seealso \code{\link{ffs}}, \code{\link{bss}}
#' @examples
#' \dontrun{
#' data(iris)
#' ffsmodel <- ffs(iris[,1:4],iris$Species)
#' plot_ffs(ffsmodel)
#' #plot performance of selected variables only:
#' plot_ffs(ffsmodel,plotType="selected")
#'}
#' @export plot_ffs
#' @aliases plot_ffs plot_bss


plot_ffs <- function(ffs_model,plotType="all",palette=rainbow,reverse=FALSE,
                     marker="black",size=1.5,lwd=0.5,
                     pch=21,...){
  message("plot_ffs() is deprecated and will be removed soon. Please use generic plot() function on ffs object.")
  plot.ffs(ffs_model, ...)


}
