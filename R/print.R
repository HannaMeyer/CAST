#' Print TrainDI
#' @param x trainDI object
#' @param ... other params
#' @export

print.trainDI = function(x, ...){
  cat(paste0("DI of ", nrow(x$train), " observation \n"))
  cat(paste0("Predictors used:"), x$variables, "\n")

  cat("\nAOA Threshold \n")
  cat(x$thres)

}


#' Show TrainDI
#' @param x trainDI object
#' @param ... other params
#' @export

show.trainDI = function(x, ...){
  print.trainDI(x)
}




#' Print AOA
#' @param x aoa object
#' @param ... other params
#' @export


print.aoa = function(x, ...){

print(x)

}


#' Show AOA
#' @param x aoa object
#' @param ... other params
#' @export


show.aoa = function(x, ...){
  print.aoa(x)
}

