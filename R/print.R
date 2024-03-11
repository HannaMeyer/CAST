#' Print CAST classes
#' @description Generic print function for trainDI and aoa
#' @name print
#' @param x trainDI object
#' @param ... other params
#' @export

print.trainDI = function(x, ...){
  cat(paste0("DI of ", nrow(x$train), " observation \n"))
  cat(paste0("Predictors:"), x$variables, "\n\n")

  cat("AOA Threshold: ")
  cat(x$threshold)

}


#' @name print
#' @param x trainDI object
#' @param ... other params
#' @export

show.trainDI = function(x, ...){
  print.trainDI(x)
}





#' @name print
#' @param x aoa object
#' @param ... other params
#' @export


print.aoa = function(x, ...){
  cat("DI:\n")
  print(x$DI)

  if ("LPD" %in% names(x)) {
    cat("LPD:\n")
    print(x$LPD)
  }

  cat("AOA:\n")
  print(x$AOA)

  cat("\n\nPredictor Weights:\n")

  print(x$parameters$weight)

  cat("\n\nAOA Threshold: ")
  cat(x$parameters$threshold)



}


#' @name print
#' @param x aoa object
#' @param ... other params
#' @export


show.aoa = function(x, ...){
  print.aoa(x)
}



#' @name print
#' @param x An object of type \emph{nndm}.
#' @param ... other arguments.
#' @export
print.nndm <- function(x, ...){
  mean_train <- round(mean(sapply(x$indx_train, length)), 2)
  min_train <- round(min(sapply(x$indx_train, length)), 2)
  cat(paste0("nndm object\n",
             "Total number of points: ", length(x$Gj), "\n",
             "Mean number of training points: ", mean_train, "\n",
             "Minimum number of training points: ", min_train, "\n"))
}

#' @name print
#' @param x An object of type \emph{nndm}.
#' @param ... other arguments.
#' @export
show.nndm = function(x, ...){
  print.nndm(x)
}


#' @name print
#' @param x An object of type \emph{knndm}.
#' @param ... other arguments.
#' @export
print.knndm <- function(x, ...){
  cat(paste0("knndm object\n",
             "Space: ", x$space, "\n",
             "Clustering algorithm: ", x$method, "\n",
             "Intermediate clusters (q): ", x$q, "\n",
             "W statistic: ", round(x$W, 4), "\n",
             "Number of folds: ", length(unique(x$clusters)),  "\n",
             "Observations in each fold: "),
      table(x$clusters), "\n")
}

#' @name print
#' @param x An object of type \emph{knndm}.
#' @param ... other arguments.
#' @export
show.knndm = function(x, ...){
  print.knndm(x)
}



#' @name print
#' @param x An object of type \emph{ffs}
#' @param ... other arguments.
#' @export


print.ffs = function(x, ...){
  cat("Selected Variables: \n")
  cat(x$selectedvars)
  cat("\n")
  cat("---\n")
  print.train(x)
}


#' @name print
#' @param x An object of type \emph{ffs}
#' @param ... other arguments.
#' @export

show.ffs = function(x, ...){
  print.ffs(x)

}




