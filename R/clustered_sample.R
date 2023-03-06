#' Clustered samples simulation
#'
#' @description A simple procedure to simulate clustered points based on a two-step sampling.
#' @param sarea polygon. Area where samples should be simulated.
#' @param nsamples integer. Number of samples to be simulated.
#' @param nparents integer. Number of parents.
#' @param radius integer. Radius of the buffer around each parent for offspring simulation.
#'
#' @return sf object with the simulated points and the parent to which each point belongs to.
#' @details A simple procedure to simulate clustered points based on a two-step sampling.
#' First, a pre-specified number of parents are simulated using random sampling.
#' For each parent, `(nsamples-nparents)/nparents` are simulated within a radius of the parent point using random sampling.
#'
#' @examples
#' # Simulate 100 points in a 100x100 square with 5 parents and a radius of 10.
#' library(sf)
#' library(ggplot2)
#'
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' simpoints <- clustered_sample(simarea, 100, 5, 10)
#' simpoints$parent <- as.factor(simpoints$parent)
#' ggplot() +
#'     geom_sf(data = simarea, alpha = 0) +
#'     geom_sf(data = simpoints, aes(col = parent))
#'
#' @author Carles Milà
#' @export
clustered_sample <- function(sarea, nsamples, nparents, radius){

  # Number of offspring per parent
  nchildren <- round((nsamples-nparents)/nparents, 0)

  # Simulate parents
  parents <- sf::st_sf(geometry=sf::st_sample(sarea, nparents, type="random"))
  res <- parents
  res$parent <- 1:nrow(parents)

  # Simulate offspring
  for(i in 1:nrow(parents)){

    # Generate buffer and cut parts outside of the area of study
    buf <- sf::st_buffer(parents[i,], dist=radius)
    buf <- sf::st_intersection(buf, sarea)

    # Simulate children
    children <- sf::st_sf(geometry=sf::st_sample(buf, nchildren, type="random"))
    children$parent <- i
    res <- rbind(res, children)
  }

  return(res)
}
