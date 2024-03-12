## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.83)

## ----message = FALSE, warning=FALSE-------------------------------------------
library(CAST)
library(caret)
library(terra)
library(sf)
library(viridis)
library(gridExtra)

## ----message = FALSE,include=FALSE, warning=FALSE-----------------------------
RMSE = function(a, b){
    sqrt(mean((a - b)^2,na.rm=T))
}

## ----message = FALSE, warning=FALSE-------------------------------------------
predictors <- rast(system.file("extdata","bioclim.tif",package="CAST"))
plot(predictors,col=viridis(100))

## ----message = FALSE, warning=FALSE-------------------------------------------

generate_random_response <- function(raster, predictornames =
names(raster), seed = sample(seq(1000), 1)){
  operands_1 = c("+", "-", "*", "/")
  operands_2 = c("^1","^2")
  
  expression <- paste(as.character(predictornames, sep=""))
  # assign random power to predictors
  set.seed(seed)
  expression <- paste(expression,
                      sample(operands_2, length(predictornames),
replace = TRUE),
                      sep = "")
  
  # assign random math function between predictors (expect after the last one)
  set.seed(seed)
  expression[-length(expression)] <- paste(expression[-
length(expression)],
                                           sample(operands_1,
length(predictornames)-1, replace = TRUE),
                                           sep = " ")
  print(paste0(expression, collapse = " "))
  # collapse
  e = paste0("raster$", expression, collapse = " ")
  
  response = eval(parse(text = e))
  names(response) <- "response"
  return(response)
  
}

## ----message = FALSE, warning=FALSE-------------------------------------------
response <- generate_random_response (predictors, seed = 10)
plot(response,col=viridis(100),main="virtual response")

## ----message = FALSE, warning=FALSE-------------------------------------------
mask <- predictors[[1]]
values(mask)[!is.na(values(mask))] <- 1
mask <- st_as_sf(as.polygons(mask))
mask <- st_make_valid(mask)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(15)
samplepoints <- st_as_sf(st_sample(mask,20,"random"))

plot(response,col=viridis(100))
plot(samplepoints,col="red",add=T,pch=3)

## ----message = FALSE, warning=FALSE-------------------------------------------
trainDat <- extract(predictors,samplepoints,na.rm=FALSE)
trainDat$response <- extract(response,samplepoints,na.rm=FALSE, ID=FALSE)$response
trainDat <- na.omit(trainDat)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model <- train(trainDat[,names(predictors)],
               trainDat$response,
               method="rf",
               importance=TRUE,
               trControl = trainControl(method="cv"))
print(model)


## ----message = FALSE, warning=FALSE-------------------------------------------
plot(varImp(model,scale = F),col="black")

## ----message = FALSE, warning=FALSE-------------------------------------------
prediction <- predict(predictors,model,na.rm=T)
truediff <- abs(prediction-response)
plot(rast(list(prediction,response)),main=c("prediction","reference"))

## ----message = FALSE, warning=FALSE-------------------------------------------
AOA <- aoa(predictors, model)
class(AOA)
names(AOA)
print(AOA)

## ----message = FALSE, warning=FALSE-------------------------------------------
plot(AOA)

## ----message = FALSE, warning=FALSE,  fig.show="hold", out.width="30%"--------
plot(truediff,col=viridis(100),main="true prediction error")
plot(AOA$DI,col=viridis(100),main="DI")
plot(prediction, col=viridis(100),main="prediction for AOA")
plot(AOA$AOA,col=c("grey","transparent"),add=T,plg=list(x="topleft",box.col="black",bty="o",title="AOA"))

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(25)
samplepoints <- clustered_sample(mask,75,15,radius=25000)

plot(response,col=viridis(100))
plot(samplepoints,col="red",add=T,pch=3)


## ----message = FALSE, warning=FALSE-------------------------------------------

trainDat <- extract(predictors,samplepoints,na.rm=FALSE)
trainDat$response <- extract(response,samplepoints,na.rm=FALSE)$response
trainDat <- data.frame(trainDat,samplepoints)
trainDat <- na.omit(trainDat)

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(10)
model_random <- train(trainDat[,names(predictors)],
               trainDat$response,
               method="rf",
               importance=TRUE,
               trControl = trainControl(method="cv"))
prediction_random <- predict(predictors,model_random,na.rm=TRUE)
print(model_random)

## ----message = FALSE, warning=FALSE-------------------------------------------
folds <- CreateSpacetimeFolds(trainDat, spacevar="parent",k=10)
set.seed(15)
model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                     method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",index=folds$index))
  print(model)
  
prediction <- predict(predictors,model,na.rm=TRUE)

## ----message = FALSE, warning=FALSE-------------------------------------------
AOA_spatial <- aoa(predictors, model)

AOA_random <- aoa(predictors, model_random)

## ----message = FALSE, warning=FALSE,  fig.show="hold", out.width="50%"--------
plot(AOA_spatial$DI,col=viridis(100),main="DI")
plot(prediction, col=viridis(100),main="prediction for AOA \n(spatial CV error applies)")
plot(AOA_spatial$AOA,col=c("grey","transparent"),add=TRUE,plg=list(x="topleft",box.col="black",bty="o",title="AOA"))
plot(prediction_random, col=viridis(100),main="prediction for AOA \n(random CV error applies)")
plot(AOA_random$AOA,col=c("grey","transparent"),add=TRUE,plg=list(x="topleft",box.col="black",bty="o",title="AOA"))

## ----message = FALSE, warning=FALSE-------------------------------------------
grid.arrange(plot(AOA_spatial) + ggplot2::ggtitle("Spatial CV"),
             plot(AOA_random) + ggplot2::ggtitle("Random CV"), ncol = 2)

## ----message = FALSE, warning=FALSE-------------------------------------------
###for the spatial CV:
RMSE(values(prediction)[values(AOA_spatial$AOA)==1],
     values(response)[values(AOA_spatial$AOA)==1])
RMSE(values(prediction)[values(AOA_spatial$AOA)==0],
     values(response)[values(AOA_spatial$AOA)==0])
model$results

###and for the random CV:
RMSE(values(prediction_random)[values(AOA_random$AOA)==1],
     values(response)[values(AOA_random$AOA)==1])
RMSE(values(prediction_random)[values(AOA_random$AOA)==0],
     values(response)[values(AOA_random$AOA)==0])
model_random$results

## ----message = FALSE, warning=FALSE-------------------------------------------
DI_RMSE_relation <- errorProfiles(model, AOA_spatial$parameters, multiCV=TRUE,
                                    window.size = 5, length.out = 5)
plot(DI_RMSE_relation)

expected_RMSE = terra::predict(AOA_spatial$DI, DI_RMSE_relation)

# account for multiCV changing the DI threshold
updated_AOA = AOA_spatial$DI > attr(DI_RMSE_relation, "AOA_threshold")


plot(expected_RMSE,col=viridis(100),main="expected RMSE")
plot(updated_AOA, col=c("grey","transparent"),add=TRUE,plg=list(x="topleft",box.col="black",bty="o",title="AOA"))

## ----message = FALSE, warning=FALSE-------------------------------------------
dat <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
# calculate average of VW for each sampling site:
dat <- aggregate(dat[,c("VW","Easting","Northing")],by=list(as.character(dat$SOURCEID)),mean)
# create sf object from the data:
pts <- st_as_sf(dat,coords=c("Easting","Northing"))

##### Extract Predictors for the locations of the sampling points
studyArea <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
st_crs(pts) <- crs(studyArea)
trainDat <- extract(studyArea,pts,na.rm=FALSE)
pts$ID <- 1:nrow(pts)
trainDat <- merge(trainDat,pts,by.x="ID",by.y="ID")
# The final training dataset with potential predictors and VW:
head(trainDat)

## ----message = FALSE, warning=FALSE-------------------------------------------
predictors <- c("DEM","NDRE.Sd","TWI","Bt")
response <- "VW"

model <- train(trainDat[,predictors],trainDat[,response],
               method="rf",tuneLength=3,importance=TRUE,
               trControl=trainControl(method="LOOCV"))
model

## ----message = FALSE, warning=FALSE-------------------------------------------
#Predictors:
plot(stretch(studyArea[[predictors]]))

#prediction:
prediction <- predict(studyArea,model,na.rm=TRUE)

## ----message = FALSE, warning=FALSE,  fig.show="hold", out.width="50%"--------
AOA <- aoa(studyArea,model)

#### Plot results:
plot(AOA$DI,col=viridis(100),main="DI with sampling locations (red)")
plot(pts,zcol="ID",col="red",add=TRUE)

plot(prediction, col=viridis(100),main="prediction for AOA \n(LOOCV error applies)")
plot(AOA$AOA,col=c("grey","transparent"),add=TRUE,plg=list(x="topleft",box.col="black",bty="o",title="AOA"))


