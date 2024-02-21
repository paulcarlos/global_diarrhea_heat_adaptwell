############## Resampling methods applied for gridded data ################
# NOTE: codes below are just the functions used, full extraction fo data not included
library(ncdf4)
library(raster)
library(sf)

# choose reference raster
ras1 <- raster("grids_country-assigned_world.tif")

# for grids smaller than 10 km; aggregation to larger grids
br1 <- raster() # put file here, raster() can be used
rs1 <- aggregate(br1,fact=10,fun=sum,expand=T) # adjust fact/factor appropriately
cellStats(rs1,"sum") # check the data loss 

# for grids greater than 10 km; resampling consistent with values within the larger grids
br1 <- raster() # put file here, raster() can be used
rs1 <- projectRaster(rs1,ras1,method="bilinear")
cellStats(rs1,"sum") # check the data loss



############## random forest for downscaling #################
library(randomForest)
library(caret)

# load data
d1 <- readRDS("combined_data_downscaling.rds")

# add to zeroes
inc.val <- 0.001

# make log rate
d1$lrate <- log(d1$drate+inc.val)

# select columns for outcome and predictors
outc <- "lrate"
#outc <- "drate" # not logged
pred <- c("yr","popden","gdppc","temp","rain","urban","crop","pasture","dist.wat","elevation")

# split data
set.seed(10272023)
ind <- sample(2,nrow(d1),replace=T,prob=c(0.7,0.3))
train <- d1[ind==1,c(outc,pred)]
test <- d1[ind==2,c(outc,pred)]
xdat <- train[,pred]
ydat <- train[,outc]

mod <- tuneRF(x=xdat,y=ydat,
              mtryStart = 3,
              ntreeTry = 1000,
              improve = 0.02,
              stepFactor = 1.5,
              trace = T,
              doBest = T, # chooses the best model, not matrix of all tried
              nodesize=length(ydat)/1000,
              importance = T,
              sampsize = length(ydat),
              replace = T,
              plot = T)
mod
pred1 <- predict(mod,test)
sqrt((mean(test$drate-pred1)^2))
sqrt((median(test$drate-pred1)^2))
plot(test$drate,pred1,xlab="actual",ylab="predicted")
varImpPlot(mod,type=1)

# save
saveRDS(mod,"model_downscale_random-forest.rds")



############## gradient boost for downscaling #################
library(caret)
library(gbm)

# load data
d1 <- readRDS("combined_data_downscaling.rds")

# add to zeroes
inc.val <- 0.001

# make log rate
d1$lrate <- log(d1$drate+inc.val)

# select columns for outcome and predictors
outc <- "lrate"
#outc <- "drate" # not logged
pred <- c("yr","popden","gdppc","temp","rain","urban","crop","pasture","dist.wat","elevation")

# split data
set.seed(10272023)
ind <- sample(2,nrow(d1),replace=T,prob=c(0.7,0.3))
train <- d1[ind==1,c(outc,pred)]
test <- d1[ind==2,c(outc,pred)]

# sample model
set.seed(10272023)
dnom <- 1e+05
train$drate <- train$drate/dnom # per hundred thousand
mod <- gbm(formula=drate~.,
           distribution="gaussian",
           data=train,
           n.trees=10000, # optimal trees minimizing the loss function
           interaction.depth=3, # splits per tree, 1 is a stump
           n.minobsinnode=15, # minimum number of observations allowed in trees
           shrinkage=0.1, # learning rate
           bag.fraction=0.65,
           train.fraction=1, # subsampling, valid.error if below 1
           cv.folds=5,
           n.cores=4,
           verbose=F)
mod
summary(mod)
which.min(mod$cv.error)
sqrt(min(mod$cv.error)) #RMSE
#sqrt(min(mod$valid.error))
#which.min(mod$valid.error)
gbm.perf(mod,method="cv")

# test predict
pred <- predict(mod,test)
v1 <- test$drate
v2 <- pred*dnom
sqrt(mean((v1-v2)^2))
sqrt(median((v1-v2)^2))
plot(v1,v2,xlab="actual",ylab="predicted")

# save
saveRDS(mod,"model_downscale_gbm.rds")



############## create proportion rasters #################
# NOTE: random forest model was used; annual rasters are in 10 year intervals

# load grid file
gr1 <- readRDS("grids_country_region_dataframe.rds")

# load model
m.rf <- readRDS("model_downscale_rf.rds")

# historical
yr <- c(2000,2010)
for (i in yr) {
  # get files
  temp <- readRDS() # temperature data from ISIMIP
  rain <- readRDS() # rainfall data from ISIMIP
  urban <- readRDS() # urban land from Hurtt et al 2020
  crop <- readRDS() # cropland from Hurtt et al 2020
  pasture <- readRDS() # pastureland from Hurtt et al 2020
  gdppc <- readRDS() # GDP from Murakami et al 2021
  popden <- readRDS() # pop density from CIESIN / Gao 2020
  pop <- readRDS() # pop counts from CIESIN / Gao 2020
  # input into a single dataframe
  gr2 <- gr1
  gr2$pop <- pop[gr2$cell]
  gr2$popden <- popden[gr2$cell]
  gr2$gdppc <- gdppc[gr2$cell]
  gr2$urban <- urban[gr2$cell]
  gr2$crop <- crop[gr2$cell]
  gr2$pasture <- pasture[gr2$cell]
  gr2$temp <- temp[gr2$cell]
  gr2$rain <- rain[gr2$cell]
  gr2$dist.wat <- dwat[gr2$cell]
  gr2$elevation <- elev[gr2$cell]
  gr2$yr <- i-2000
  #sum(is.na(gr2))
  gr2[is.na(gr2)] <- 0
  # run models 
  gr2$drate.rf <- predict(m.rf,gr2)
  # save file
  saveRDS(gr2,paste0("grids_",i,"_dataframe.rds"))
  # remove files to increase RAM
  rm(temp,rain,urban,crop,pasture,gdppc,popden,pop,gr2);gc()
}
rm(i)

# future
ssp <- c("ssp126","ssp370","ssp585")
yr <- seq(2020,2100,10)
for (i in seq(ssp)) {
  s1 <- ssp[i]
  for (j in seq(yr)) {
    y1 <- yr[j]
    # get files
    temp <- readRDS() # temperature data from ISIMIP
    rain <- readRDS() # rainfall data from ISIMIP
    urban <- readRDS() # urban land from Hurtt et al 2020
    crop <- readRDS() # cropland from Hurtt et al 2020
    pasture <- readRDS() # pastureland from Hurtt et al 2020
    gdppc <- readRDS() # GDP from Murakami et al 2021
    popden <- readRDS() # pop density from CIESIN / Gao 2020
    pop <- readRDS() # pop counts from CIESIN / Gao 2020
    # input into a single dataframe
    gr2 <- gr1
    gr2$pop <- pop[gr2$cell]
    gr2$popden <- gr2$pop/85 # km^2 assuming 0.083 degrees is 9.2 km
    gr2$gdppc <- gdp[gr2$cell]/pop[gr2$cell]
    gr2$urban <- urban[gr2$cell]
    gr2$crop <- crop[gr2$cell]
    gr2$pasture <- pasture[gr2$cell]
    gr2$temp <- temp[gr2$cell]
    gr2$rain <- rain[gr2$cell]
    gr2$dist.wat <- dwat[gr2$cell]
    gr2$elevation <- elev[gr2$cell]
    gr2$yr <- y1-2000
    #sum(is.na(gr2))
    gr2[is.na(gr2)] <- 0
    gr2$gdppc[is.infinite(gr2$gdppc)] <- 0
    #summary(gr2)
    # run models 
    gr2$drate.rf <- predict(m.rf,gr2)
    # save file
    saveRDS(gr2,paste0("grids_",s1,"_",y1,"_dataframe.rds"))
    # remove files to increase RAM
    rm(temp,rain,urban,crop,pasture,gdp,pop,gr2);gc()
    cat(paste0(i,"-",j)," ")
  }
}
rm(i,j,s1,y1)
