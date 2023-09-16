#' [EXPERIMENTAL] K-fold Nearest Neighbour Distance Matching
#' @description
#' This function implements the kNNDM algorithm and returns the necessary
#' indices to perform a k-fold NNDM CV for map validation.
#'
#' @author Carles Milà and Jan Linnenbrink
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param modeldomain sf polygon object defining the prediction area. Optional; alternative to ppoints (see Details).
#' @param ppoints sf or sfc point object. Contains the target prediction points. Optional; alternative to modeldomain (see Details).
#' @param space character. Only "geographical" knndm, i.e. kNNDM in the geographical space, is currently implemented.
#' @param k integer. Number of folds desired for CV. Defaults to 10.
#' @param maxp numeric. Maximum fold size allowed, defaults to 0.5, i.e. a single fold can hold a maximum of half of the training points.
#' @param clustering character. Possible values include "hierarchical" and "kmeans". See details.
#' @param linkf character. Only relevant if clustering = "hierarchical". Link function for agglomerative hierarchical clustering.
#' Defaults to "ward.D2". Check `stats::hclust` for other options.
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of ppoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See `sf::st_sample`.
#' Only required if modeldomain is used instead of ppoints.
#'
#' @return An object of class \emph{knndm} consisting of a list of eight elements:
#' indx_train, indx_test (indices of the observations to use as
#' training/test data in each kNNDM CV iteration), Gij (distances for
#' G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during kNNDM CV), clusters (list of cluster IDs),
#' W (Wasserstein statistic), and space (stated by the user in the function call).
#'
#' @details
#' knndm is a k-fold version of NNDM LOO CV for medium and large datasets. Brielfy, the algorithm tries to
#' find a k-fold configuration such that the integral of the absolute differences (Wasserstein W statistic)
#' between the empirical nearest neighbour distance distribution function between the test and training data during CV (Gj*),
#' and the empirical nearest neighbour distance distribution function between the prediction and training points (Gij),
#' is minimised. It does so by performing clustering of the training points' coordinates for different numbers of
#' clusters that range from k to N (number of observations), merging them into k final folds,
#' and selecting the configuration with the lowest W.
#'
#' Using a projected CRS in `knndm` has large computational advantages since fast nearest neighbour search can be
#' done via the `FNN` package, while working with geographic coordinates requires computing the full
#' spherical distance matrices. As a clustering algorithm, `kmeans` can only be used for
#' projected CRS while `hierarchical` can work with both projected and geographical coordinates, though it requires
#' calculating the full distance matrix of the training points even for a projected CRS.
#'
#' In order to select between clustering algorithms and number of folds `k`, different `knndm` configurations can be run
#' and compared, being the one with a lower W statistic the one that offers a better match. W statistics between `knndm`
#' runs are comparable as long as `tpoints` and `ppoints` or `modeldomain` stay the same.
#'
#' Map validation using knndm should be used using `CAST::global_validation`, i.e. by stacking all out-of-sample
#' predictions and evaluating them all at once. The reasons behind this are 1) The resulting folds can be
#' unbalanced and 2) nearest neighbour functions are constructed and matched using all CV folds simultaneously.
#'
#' If training data points are very clustered with respect to the prediction area and the presented knndm
#' configuration still show signs of Gj* > Gij, there are several things that can be tried. First, increase
#' the `maxp` parameter; this may help to control for strong clustering (at the cost of having unbalanced folds).
#' Secondly, decrease the number of final folds `k`, which may help to have larger clusters.
#'
#' The `modeldomain` is a sf polygon that defines the prediction area. The function takes a regular point sample
#' (amount defined by `samplesize`) from the spatial extent. As an alternative use `ppoints` instead of
#' `modeldomain`, if you have already defined the prediction locations (e.g. raster pixel centroids).
#' When using either `modeldomain` or `ppoints`, we advise to plot the study area polygon and the
#' training/prediction points as a previous step to ensure they are aligned.
#'
#' @note Experimental cycle. Article describing and testing the algorithm in preparation.
#' @references
#' \itemize{
#' \item Linnenbrink, J., Milà, C., Ludwig, M., and Meyer, H.: kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1308, 2023.
#' \item Milà, C., Mateu, J., Pebesma, E., Meyer, H. (2022): Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation. Methods in Ecology and Evolution 00, 1– 13.
#' }
#' @seealso \code{\link{geodist}}, \code{\link{nndm}}
#'
#' @export
#' @examples
#' ########################################################################
#' # Example 1: Simulated data - Randomly-distributed training points
#' ########################################################################
#'
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 1000 random training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- sf::st_sample(simarea, 1000, type = "random")
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#' plot(simarea)
#' plot(pred_points, add = TRUE, col = "blue")
#' plot(train_points, add = TRUE, col = "red")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_folds <- knndm(train_points, ppoints = pred_points, k = 5)
#' knndm_folds
#' plot(knndm_folds)
#' folds <- as.character(knndm_folds$clusters)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, aes(col = folds))
#'
#' ########################################################################
#' # Example 2: Simulated data - Clustered training points
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(ggplot2)
#'
#' # Simulate 1000 clustered training points in a 100x100 square
#' set.seed(1234)
#' simarea <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' simarea <- sf::st_polygon(simarea)
#' train_points <- clustered_sample(simarea, 1000, 50, 5)
#' pred_points <- sf::st_sample(simarea, 1000, type = "regular")
#' plot(simarea)
#' plot(pred_points, add = TRUE, col = "blue")
#' plot(train_points, add = TRUE, col = "red")
#'
#' # Run kNNDM for the whole domain, here the prediction points are known.
#' knndm_folds <- knndm(train_points, ppoints = pred_points, k = 5)
#' knndm_folds
#' plot(knndm_folds)
#' folds <- as.character(knndm_folds$clusters)
#' ggplot() +
#'   geom_sf(data = simarea, alpha = 0) +
#'   geom_sf(data = train_points, aes(col = folds))
#'}
#' ########################################################################
#' # Example 3: Real- world example; using a modeldomain instead of previously
#' # sampled prediction locations
#' ########################################################################
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#'
#' ### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
#'    by=list(as.character(dat$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' studyArea <- rast(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#' studyArea <- as.polygons(studyArea, values = FALSE, na.all = TRUE) |>
#'     st_as_sf() |>
#'     st_union()
#' pts <- st_transform(pts, crs = st_crs(studyArea))
#' plot(studyArea)
#' plot(st_geometry(pts), add = TRUE, col = "red")
#'
#' knndm_folds <- knndm(pts, modeldomain=studyArea, k = 5)
#' knndm_folds
#' plot(knndm_folds)
#' folds <- as.character(knndm_folds$clusters)
#' ggplot() +
#'   geom_sf(data = pts, aes(col = folds))
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=knndm_folds$indx_train,
#'    savePredictions='final')
#' model_knndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
#'    dat$VW,
#'    method="rf",
#'    trControl = ctrl)
#' global_validation(model_knndm)
#'}
knndm <- function(tpoints, modeldomain = NULL, ppoints = NULL,
                  space = "geographical",
                  k = 10, maxp = 0.5,
                  clustering = "hierarchical", linkf = "ward.D2",
                  samplesize = 1000, sampling = "regular"){

  # create sample points from modeldomain
  if(is.null(ppoints)&!is.null(modeldomain)){
    if(!identical(sf::st_crs(tpoints), sf::st_crs(modeldomain))){
      stop("tpoints and modeldomain must have the same CRS")
    }
    message(paste0(samplesize, " prediction points are sampled from the modeldomain"))
    ppoints <- sf::st_sample(x = modeldomain, size = samplesize, type = sampling)
    sf::st_crs(ppoints) <- sf::st_crs(modeldomain)
  }

  # Conditional preprocessing actions
  if (any(class(tpoints) %in% "sfc")) {
    tpoints <- sf::st_sf(geom = tpoints)
  }
  if (any(class(ppoints) %in% "sfc")) {
    ppoints <- sf::st_sf(geom = ppoints)
  }
  if(is.na(sf::st_crs(tpoints))){
    warning("Missing CRS in training or prediction points. Assuming projected CRS.")
    islonglat <- FALSE
  }else{
    islonglat <- sf::st_is_longlat(tpoints)
  }

  # Prior checks
  check_knndm(tpoints, ppoints, space, k, maxp, clustering, islonglat)

  # kNNDM in the geographical space (currently only option)
  if(isTRUE(space == "geographical")){
    knndm_res <- knndm_geo(tpoints, ppoints, k, maxp, clustering, linkf, islonglat)
  }

  # Output
  knndm_res
}



# kNNDM checks
check_knndm <- function(tpoints, ppoints, space, k, maxp, clustering, islonglat){

  if(!identical(sf::st_crs(tpoints), sf::st_crs(ppoints))){
    stop("tpoints and ppoints must have the same CRS")
  }
  if (!(clustering %in% c("kmeans", "hierarchical"))) {
    stop("clustering must be one of `kmeans` or `hierarchical`")
  }
  if (space != "geographical") {
    stop("Only kNNDM in the geographical space is currently implemented.")
  }
  if (!(maxp < 1 & maxp > 1/k)) {
    stop("maxp must be strictly between 1/k and 1")
  }
  if(isTRUE(islonglat) & clustering == "kmeans"){
    stop("kmeans works in the Euclidean space and therefore can only handle
         projected coordinates. Please use hierarchical clustering or project your data.")
  }
}

# kNNDM in the geographical space
knndm_geo <- function(tpoints, ppoints, k, maxp, clustering, linkf, islonglat){

  # Gj and Gij calculation
  tcoords <- sf::st_coordinates(tpoints)[,1:2]
  if(isTRUE(islonglat)){
    distmat <- sf::st_distance(tpoints)
    units(distmat) <- NULL
    diag(distmat) <- NA
    Gj <- apply(distmat, 1, function(x) min(x, na.rm=TRUE))
    Gij <- sf::st_distance(ppoints, tpoints)
    units(Gij) <- NULL
    Gij <- apply(Gij, 1, min)
  }else{
    Gj <- c(FNN::knn.dist(tcoords, k = 1))
    Gij <- c(FNN::knnx.dist(query = sf::st_coordinates(ppoints)[,1:2],
                            data = tcoords, k = 1))
  }

  # Check if Gj > Gij (warning suppressed regarding ties)
  testks <- suppressWarnings(stats::ks.test(Gj, Gij, alternative = "great"))
  if(testks$p.value >= 0.05){

    clust <- sample(rep(1:k, ceiling(nrow(tpoints)/k)), size = nrow(tpoints), replace=F)

    if(isTRUE(islonglat)){
      Gjstar <- distclust_geo(distmat, clust)
    }else{
      Gjstar <- distclust_proj(tcoords, clust)
    }
    k_final <- "random CV"
    W_final <- twosamples::wass_stat(Gjstar, Gij)
    message("Gij <= Gj; a random CV assignment is returned")

  }else{

    if(clustering == "hierarchical"){
      # For hierarchical clustering we need to compute the full distance matrix,
      # but we can integrate geographical distances
      if(!isTRUE(islonglat)){
        distmat <- sf::st_distance(tpoints)
      }
      hc <- stats::hclust(d = stats::as.dist(distmat), method = linkf)
    }

    # Build grid of number of clusters to try - we sample low numbers more intensively
    clustgrid <- data.frame(nk = as.integer(round(exp(seq(log(k), log(nrow(tpoints)-2),
                                                          length.out = 100)))))
    clustgrid$W <- NA
    clustgrid <- clustgrid[!duplicated(clustgrid$nk),]
    clustgroups <- list()

    # Compute 1st PC for ordering clusters
    pcacoords <- stats::prcomp(tcoords, center = TRUE, scale. = FALSE, rank = 1)

    # We test each number of clusters
    for(nk in clustgrid$nk){

      # Create nk clusters
      if(clustering == "hierarchical"){
        clust_nk <- stats::cutree(hc, k=nk)
      }else if(clustering == "kmeans"){
        clust_nk <- stats::kmeans(tcoords, nk)$cluster
      }

      tabclust <- as.data.frame(table(clust_nk))
      tabclust$clust_k <- NA

      # compute cluster centroids and apply PC loadings to shuffle along the 1st dimension
      centr_tpoints <- sapply(tabclust$clust_nk, function(x){
        centrpca <- matrix(apply(tcoords[clust_nk %in% x, , drop=FALSE], 2, mean), nrow = 1)
        colnames(centrpca) <- colnames(tcoords)
        return(predict(pcacoords, centrpca))
      })
      tabclust$centrpca <- centr_tpoints
      tabclust <- tabclust[order(tabclust$centrpca),]

      # We don't merge big clusters
      clust_i <- 1
      for(i in 1:nrow(tabclust)){
        if(tabclust$Freq[i] >= nrow(tpoints)/k){
          tabclust$clust_k[i] <- clust_i
          clust_i <- clust_i + 1
        }
      }
      rm("clust_i")

      # And we merge the remaining into k groups
      clust_i <- setdiff(1:k, unique(tabclust$clust_k))
      tabclust$clust_k[is.na(tabclust$clust_k)] <- rep(clust_i, ceiling(nk/length(clust_i)))[1:sum(is.na(tabclust$clust_k))]
      tabclust2 <- data.frame(ID = 1:length(clust_nk), clust_nk = clust_nk)
      tabclust2 <- merge(tabclust2, tabclust, by = "clust_nk")
      tabclust2 <- tabclust2[order(tabclust2$ID),]
      clust_k <- tabclust2$clust_k

      # Compute W statistic if not exceeding maxp
      if(!any(table(clust_k)/length(clust_k)>maxp)){

        if(isTRUE(islonglat)){
          Gjstar_i <- distclust_geo(distmat, clust_k)
        }else{
          Gjstar_i <- distclust_proj(tcoords, clust_k)
        }
        clustgrid$W[clustgrid$nk==nk] <- twosamples::wass_stat(Gjstar_i, Gij)
        clustgroups[[paste0("nk", nk)]] <- clust_k
      }
    }

    # Final configuration
    k_final <- clustgrid$nk[which.min(clustgrid$W)]
    W_final <- min(clustgrid$W, na.rm=T)
    clust <- clustgroups[[paste0("nk", k_final)]]
    if(isTRUE(islonglat)){
      Gjstar <- distclust_geo(distmat, clust)
    }else{
      Gjstar <- distclust_proj(tcoords, clust)
    }
  }

  # Output
  cfolds <- CAST::CreateSpacetimeFolds(data.frame(clust=clust), spacevar = "clust", k = k)
  res <- list(clusters = clust,
              indx_train = cfolds$index, indx_test = cfolds$indexOut,
              Gij = Gij, Gj = Gj, Gjstar = Gjstar,
              W = W_final, method = clustering, q = k_final, space = "geographical")
  class(res) <- c("knndm", "list")
  res
}

# Helper function: Compute out-of-fold NN distance (geographical)
distclust_geo <- function(distm, folds){
  alldist <- rep(NA, length(folds))
  for(f in unique(folds)){
    alldist[f == folds] <- apply(distm[f == folds, f != folds, drop=FALSE], 1, min)
  }
  alldist
}

# Helper function: Compute out-of-fold NN distance (projected)
distclust_proj <- function(tr_coords, folds){
  alldist <- rep(NA, length(folds))
  for(f in unique(folds)){
    alldist[f == folds] <- c(FNN::knnx.dist(query = tr_coords[f == folds,,drop=FALSE],
                                            data = tr_coords[f != folds,,drop=FALSE], k = 1))
  }
  alldist
}
