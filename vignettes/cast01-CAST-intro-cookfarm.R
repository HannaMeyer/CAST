## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.83,cache = FALSE)
user_hanna <- Sys.getenv("USER") %in% c("hanna")

## ----c1, message = FALSE, warning=FALSE---------------------------------------
#install.packages("CAST")
library(CAST)

## ----c2, message = FALSE, warning=FALSE---------------------------------------
help(CAST)

## ----c3, message = FALSE, warning=FALSE---------------------------------------
data <- readRDS(system.file("extdata","Cookfarm.RDS",package="CAST"))
head(data)

## ----c4, message = FALSE, warning=FALSE---------------------------------------

library(sf)
data_sp <- unique(data[,c("SOURCEID","Easting","Northing")])
data_sp <- st_as_sf(data_sp,coords=c("Easting","Northing"),crs=26911)
plot(data_sp,axes=T,col="black")

## ----c5, message = FALSE, warning=FALSE, eval=user_hanna----------------------
#...or plot the data with mapview:
library(mapview)
mapviewOptions(basemaps = c("Esri.WorldImagery"))
mapview(data_sp)

## ----c6, message = FALSE, warning=FALSE---------------------------------------
library(lubridate)
library(ggplot2)
trainDat <- data[data$altitude==-0.3&
                   year(data$Date)==2012&
                   week(data$Date)%in%c(10:12),]
ggplot(data = trainDat, aes(x=Date, y=VW)) +
  geom_line(aes(colour=SOURCEID))

## ----c7, message = FALSE, warning=FALSE---------------------------------------
library(caret)
predictors <- c("DEM","TWI","Precip_cum","cday",
                "MaxT_wrcc","Precip_wrcc","BLD",
                "Northing","Easting","NDRE.M")
set.seed(10)
model <- train(trainDat[,predictors],trainDat$VW,
               method="rf",tuneGrid=data.frame("mtry"=2),
               importance=TRUE,ntree=50,
               trControl=trainControl(method="cv",number=3))

## ----c8, message = FALSE, warning=FALSE---------------------------------------
library(terra)
predictors_sp <- rast(system.file("extdata","predictors_2012-03-25.tif",package="CAST"))
prediction <- predict(predictors_sp,model,na.rm=TRUE)
plot(prediction)

## ----c9, message = FALSE, warning=FALSE---------------------------------------
model

## ----c10, message = FALSE, warning=FALSE--------------------------------------
set.seed(10)
indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID",
                                k=3)
set.seed(10)
model_LLO <- train(trainDat[,predictors],trainDat$VW,
                   method="rf",tuneGrid=data.frame("mtry"=2), importance=TRUE,
                   trControl=trainControl(method="cv",
                                          index = indices$index))
model_LLO

## ----c11, message = FALSE, warning=FALSE--------------------------------------
plot(varImp(model_LLO))

## ----c12, message = FALSE, warning=FALSE--------------------------------------
set.seed(10)
ffsmodel_LLO <- ffs(trainDat[,predictors],trainDat$VW,metric="Rsquared",
                    method="rf", tuneGrid=data.frame("mtry"=2),
                    verbose=FALSE,ntree=50,
                    trControl=trainControl(method="cv",
                                           index = indices$index))
ffsmodel_LLO
ffsmodel_LLO$selectedvars

## ----c13, message = FALSE, warning=FALSE--------------------------------------
plot(ffsmodel_LLO)

## ----c14, message = FALSE, warning=FALSE--------------------------------------
prediction_ffs <- predict(predictors_sp,ffsmodel_LLO,na.rm=TRUE)
plot(prediction_ffs)

## ----c15, message = FALSE, warning=FALSE--------------------------------------
### AOA for which the spatial CV error applies:
AOA <- aoa(predictors_sp,ffsmodel_LLO)

plot(prediction_ffs,main="prediction for the AOA \n(spatial CV error applied)")
plot(AOA$AOA,col=c("grey","transparent"),add=T)

#spplot(prediction_ffs,main="prediction for the AOA \n(spatial CV error applied)")+
#spplot(AOA$AOA,col.regions=c("grey","transparent"))

### AOA for which the random CV error applies:
AOA_random <- aoa(predictors_sp,model)
plot(prediction,main="prediction for the AOA \n(random CV error applied)")
plot(AOA_random$AOA,col=c("grey","transparent"),add=T)

#spplot(prediction,main="prediction for the AOA \n(random CV error applied)")+
#spplot(AOA_random$AOA,col.regions=c("grey","transparent"))


