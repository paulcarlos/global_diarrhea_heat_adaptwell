################ create baseline temperature data ##################
# NOTE: baseline years aligned with CMIP6; data from ISIMIP3b
library(raster)
library(ncdf4)

# years
yr <- 1985:2014 # 30 years

# model
#m1 <- "ipsl-cm6a-lr"
m1 <- "mri-esm2-0"

# files
fn1 <- list.files() # get files

# create full array
r1 <- readRDS(fn1[1])
for (i in 2:length(fn1)) {
  r2 <- readRDS(fn1[i])
  r1 <- addLayer(r1,r2)
}
rm(i,r2)

r3 <- calc(r1,fun=mean,na.rm=T)
plot(r3)

# save
saveRDS(r3,paste0("temp_reference_",m1,"_1985-2014.rds"))


################ create temperature difference data ##################
library(raster)
#library(weathermetrics)

ras1 <- raster("") # for resampling
f1 <- "" # folder containing annual files

# model
#m1 <- "ipsl-cm6a-lr"
m1 <- "mri-esm2-0"

# get reference files
rf <- readRDS(paste0("reference/temp_reference_",m1,"_1985-2014.rds"))

# loop historical
yr <- 2000:2014
for (i in yr) {
  d1 <- readRDS(list.files(paste0(f1,"annual/",m1,"/historical/"),pattern=paste0("_",i),full.names=T))
  d2 <- d1-rf
  d3 <- projectRaster(d2,ras1,method="ngb")
  writeRaster(d3,paste0(f1,"difference/",m1,"/historical/temp_diff_",m1,"_",i,".tif"),overwrite=T)
  rm(d1,d2,d3);gc()
  cat(i," ")
}
rm(i,yr)
yr <- 2015:2019 #use SSP2
for (i in yr) {
  d1 <- readRDS(list.files(paste0(f1,"annual/",m1,"/projected/"),pattern=paste0("ssp245_",i),full.names=T))
  d2 <- d1-rf
  d3 <- projectRaster(d2,ras1,method="ngb")
  writeRaster(d3,paste0(f1,"difference/",m1,"/historical/temp_diff_",m1,"_",i,".tif"),overwrite=T)
  rm(d1,d2,d3);gc()
  cat(i," ")
}
rm(i,yr)

# loop projected
ssp <- c("ssp126","ssp370","ssp585")
yr <- 2020:2100
for (i in ssp) {
  for (j in yr) {
    d1 <- readRDS(list.files(paste0(f1,"annual/",m1,"/projected/"),pattern=paste0(i,"_",j),full.names=T))
    d2 <- d1-rf
    d3 <- projectRaster(d2,ras1,method="ngb")
    writeRaster(d3,paste0(f1,"difference/",m1,"/projected/temp_diff_",m1,"_",i,"_",j,".tif"),overwrite=T)
    rm(d1,d2,d3);gc()
    cat(paste0(i,"-",j)," ") 
  }
}
rm(i,j)
