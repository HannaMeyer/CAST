## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.width=6.2, fig.height=3.4)


## ----message = FALSE, warning=FALSE-------------------------------------------
library(CAST)
library(caret)
library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)

## ----message = FALSE, warning=FALSE-------------------------------------------
seed <- 10 # random realization
samplesize <- 300 # how many samples will be used?
nparents <- 20 #For clustered samples: How many clusters? 
radius <- 500000 # For clustered samples: What is the radius of a cluster?


## ----message = FALSE, warning=FALSE-------------------------------------------
ee <- st_crs("+proj=eqearth")
co <- ne_countries(returnclass = "sf")
co.ee <- st_transform(co, ee)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
sf_use_s2(FALSE)
set.seed(seed)
pts_random <- st_sample(co.ee, samplesize)
### See points on the map:
ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_random, color = "#F8766D",size=0.5, shape=3) +
  guides(fill = "none", col = "none") +
  labs(x = NULL, y = NULL)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
set.seed(seed)
sf_use_s2(FALSE)
pts_clustered <- clustered_sample(co.ee, samplesize, nparents, radius)

ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = "#F8766D",size=0.5, shape=3) +
  guides(fill = "none", col = "none") +
  labs(x = NULL, y = NULL)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_random <- geodist(pts_random,co.ee,
                            sampling="Fibonacci")
dist_clstr <- geodist(pts_clustered,co.ee,
                           sampling="Fibonacci")

plot(dist_random, unit = "km")+scale_x_log10(labels=round)+ggtitle("Randomly distributed reference data")
plot(dist_clstr, unit = "km")+scale_x_log10(labels=round)+ggtitle("Clustered reference data")
             

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
randomfolds <- caret::createFolds(1:nrow(pts_clustered))

## ----message = FALSE, warning=FALSE, results='hide',echo=FALSE----------------
for (i in 1:nrow(pts_clustered)){
  pts_clustered$randomCV[i] <- which(unlist(lapply(randomfolds,function(x){sum(x%in%i)}))==1)
}

ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = rainbow(max(pts_clustered$randomCV))[pts_clustered$randomCV],size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ggtitle("random fold membership shown by color")

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_clstr <- geodist(pts_clustered,co.ee,
                           sampling="Fibonacci", 
                           cvfolds= randomfolds)
plot(dist_clstr, unit = "km")+scale_x_log10(labels=round)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
spatialfolds <- CreateSpacetimeFolds(pts_clustered,spacevar="parent",k=length(unique(pts_clustered$parent)))

## ----message = FALSE, warning=FALSE, results='hide',echo=FALSE----------------
ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, color = rainbow(max(pts_clustered$parent))[pts_clustered$parent],size=0.5, shape=3) +
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ ggtitle("spatial fold membership by color")

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
dist_clstr <- geodist(pts_clustered,co.ee,
                           sampling="Fibonacci",
                           cvfolds= spatialfolds$indexOut)
plot(dist_clstr, unit = "km")+scale_x_log10(labels=round)

             

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
# create a spatial CV for the randomly distributed data. Here:
# "leave region-out-CV"
sf_use_s2(FALSE)
pts_random_co <- st_join(st_as_sf(pts_random),co.ee)


ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_random_co, aes(color=subregion),size=0.5, shape=3) +
  scale_color_manual(values=rainbow(length(unique(pts_random_co$subregion))))+
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ ggtitle("spatial fold membership by color")


## ----message = FALSE, warning=FALSE, results='hide'---------------------------

spfolds_rand <- CreateSpacetimeFolds(pts_random_co,spacevar = "subregion",
                                     k=length(unique(pts_random_co$subregion)))
dist_rand_sp <- geodist(pts_random_co,co.ee,
                             sampling="Fibonacci", 
                             cvfolds= spfolds_rand$indexOut)
plot(dist_rand_sp, unit = "km")+scale_x_log10(labels=round)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------


nndmfolds_clstr <- nndm(pts_clustered, modeldomain=co.ee, samplesize = 2000)
dist_clstr <- geodist(pts_clustered,co.ee,
                           sampling = "Fibonacci",
                           cvfolds = nndmfolds_clstr$indx_test, 
                           cvtrain = nndmfolds_clstr$indx_train)
plot(dist_clstr, unit = "km")+scale_x_log10(labels=round)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------

nndmfolds_rand <- nndm(pts_random_co,  modeldomain=co.ee, samplesize = 2000)
dist_rand <- geodist(pts_random_co,co.ee,
                          sampling = "Fibonacci",
                          cvfolds = nndmfolds_rand$indx_test, 
                          cvtrain = nndmfolds_rand$indx_train)
plot(dist_rand, unit = "km")+scale_x_log10(labels=round)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------

knndmfolds_clstr <- knndm(pts_clustered, modeldomain=co.ee, samplesize = 2000)
pts_clustered$knndmCV <- as.character(knndmfolds_clstr$clusters)

ggplot() + geom_sf(data = co.ee, fill="#00BFC4",col="#00BFC4") +
  geom_sf(data = pts_clustered, aes(color=knndmCV),size=0.5, shape=3) +
  scale_color_manual(values=rainbow(length(unique(pts_clustered$knndmCV))))+
  guides(fill = FALSE, col = FALSE) +
  labs(x = NULL, y = NULL)+ ggtitle("spatial fold membership by color")


dist_clstr <- geodist(pts_clustered,co.ee,
                           sampling = "Fibonacci",
                           cvfolds = knndmfolds_clstr$indx_test, 
                           cvtrain = knndmfolds_clstr$indx_train)
plot(dist_clstr, unit = "km")+scale_x_log10(labels=round)


## ----message = FALSE, warning=FALSE, results='hide'---------------------------
predictors_global <- rast(system.file("extdata","bioclim_global.tif",package="CAST"))

plot(predictors_global)

## ----message = FALSE, warning=FALSE, results='hide'---------------------------

# use random CV:
dist_clstr_rCV <- geodist(pts_clustered,predictors_global,
                               type = "feature", 
                               sampling="Fibonacci",
                               cvfolds = randomfolds)

# use spatial CV:
dist_clstr_sCV <- geodist(pts_clustered,predictors_global,
                               type = "feature", sampling="Fibonacci",
                               cvfolds = spatialfolds$indexOut)


# Plot results:
plot(dist_clstr_rCV)+scale_x_log10()+ggtitle("Clustered reference data and random CV")
plot(dist_clstr_sCV)+scale_x_log10()+ggtitle("Clustered reference data and spatial CV")

