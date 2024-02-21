################### calculating excess diarrheal deaths ####################
# NOTE: for same adaptation in 2019 in the future; codes are similar with other adaptation options

library(raster)
library(lme4)
library(ncdf4)

# load countries 
ctry1 <- read.csv("locations_compiled_ihme_wcde_pik_unicef.csv")
ctry2 <- ctry1[!is.na(ctry1$wcde2.ctry) & ctry1$ihme.iso3!="TWN",] # updated data from PIK removed Taiwan
loc <- sort(ctry2$ihme.iso3)

# load data
ihme.mod <- readRDS("lmer_model_w-paf_6-scalar.rds") # model
sdi <- readRDS("sdi_1990-2100_ssp1-5.rds") # SDI data
pop <- readRDS("wcde_pop_hist-proj_1990-2100_ssp1-5_ihme-countries.rds") # population in thousands
pop$age <- gsub("--","-",pop$age)
gr <- readRDS("grids_country_region_dataframe.rds") # grids
scal <- readRDS("risk_factors_log-scalar_2019.rds")
ras <- raster("grids_country-assigned_world.tif")
ras[] <- NA # empty
rr <- read.csv("pct-change_pathogen_chua2021.csv")

# file locations
fsav <- "" # folder to save files
fras <- "" # folder for downscaled proportion rasters
ftmp <- "" # folder for temperatures

# categories
yr1 <- c("2020-2039","2040-2059","2060-2079","2080-2099")
yr2 <- list(2020:2039,2040:2059,2060:2079,2080:2099)
yref <- 2000
ssp.reg <- c("SSP1","SSP3","SSP5")
ssp.tmp <- c("ssp126","ssp370","ssp585")
rei1 <- c("camp","chol","cryp","epec","etec","noro","rota","salm","shig","typh")
rei2 <- c("Campylobacter","Cholera","Cryptosporidium","Enteropathogenic E coli","Enterotoxigenic E coli",
          "Norovirus","Rotavirus","Non-typhoidal Salmonella","Shigella","Typhoid fever")
age <- c("0-4","5-69","70+")
#acat <- c("under5","5plus")
inc.val <- 0.001



######### loop codes
r1 <- "mean" # options are: "mean" "lower" "upper"
nm2 <- "scal.o5.all" # for >5 yo
dg <- 6 # digits to round off
#i=1;j=2;k=4
for (i in seq(rei1)) {
  p1 <- rei1[i] # select pathogen
  p2 <- rei2[i]
  if (p1=="rota") {
    nm1 <- "scal.u5.rota"
  } else {
    nm1 <- "scal.u5.all"
  }
  for (j in seq(ssp.reg)) {
    s1 <- ssp.reg[j] # SSP scenarios
    s2 <- ssp.tmp[j]
    for (k in seq(yr1)) {
      y1 <- yr1[k] # year category
      y2 <- yr2[[k]]
      
      # choose rasters
      pr.u5 <- pr.o5 <- ras # blank raster for proportion
      ras.pr <- raster(paste0(fras,"grids_",tolower(s1),"_",y1,"_prop_raster.tif"))
      tmp1 <- raster(paste0(ftmp,"temp_diff_model-av_",s2,"_",y1,".tif"))
      rval <- rr[rr$rei==p1,r1]
      tmp2 <- exp(tmp1*rval)
      temp <- (tmp2-1)/tmp2
      #cellStats(temp,"range")
      rm(tmp1,tmp2);gc()
      
      # get population and SDI data
      sdi1 <- sdi[sdi$ssp==s1,]
      pop1 <- pop[pop$ssp==s1,]
      
      # create base dataframe
      df1 <- expand.grid("country"=loc,"age"=age,"year"=y2)
      df1$yr <- df1$year-yref
      df1$pathogen <- p2
      df1$sdi <- sdi1$sdi[match(paste0(df1$country,"-",df1$year),paste0(sdi1$iso3,"-",sdi1$year))]
      df1$pop <- pop1$pop[match(paste0(df1$country,"-",df1$age,"-",df1$year),
                                paste0(pop1$iso3,"-",pop1$age,"-",pop1$year))]
      df1$logscal <- NA
      df1$logscal[df1$age=="0-4"] <- scal[match(df1$country[df1$age=="0-4"],scal$iso3),nm1]
      df1$logscal[df1$age!="0-4"] <- scal[match(df1$country[df1$age!="0-4"],scal$iso3),nm2]
      
      # get counts
      df1$drate <- exp(predict(ihme.mod,newdata=df1))-inc.val
      df1$counts <- (df1$drate/100)*df1$pop # per 100 thousand
      
      # get average
      df.u5 <- aggregate(counts~country,df1[df1$age=="0-4",],"mean")
      df.o5 <- aggregate(counts~country+year,df1[df1$age!="0-4",],"sum")
      df.o5 <- aggregate(counts~country,df.o5,"mean")
      #sum(df.u5$counts)
      #sum(df.o5$counts)
      
      # put in grids
      gr1 <- gr[,c("cell","iso3")]
      gr1$u5 <- df.u5$counts[match(gr1$iso3,df.u5$country)]
      gr1$o5 <- df.o5$counts[match(gr1$iso3,df.o5$country)]
      #gr1$prop <- ras.pr[gr1$cell]
      pr.u5[gr1$cell] <- gr1$u5
      pr.o5[gr1$cell] <- gr1$o5
      #ras.u5 <- pr.u5*ras.pr*temp
      #ras.o5 <- pr.o5*ras.pr*temp
      ras.u5 <- round(pr.u5*ras.pr*temp,dg)
      ras.o5 <- round(pr.o5*ras.pr*temp,dg)
      #plot(ras.u5)
      #plot(ras.o5)
      
      #save
      writeRaster(ras.u5, 
                  paste0(fsav,"temp-deaths_under5_",p1,"_",s2,"_",y1,"_",r1,"_baseline2019.nc"), #filename 
                  overwrite=T, format="CDF", datatype="FLT4S", force_v4=TRUE, compression=5,  
                  varname=paste0("mean-temp-mort_u5_",p1), # variable name
                  varunit="human deaths", 
                  longname=paste0("mean temperature-related excess deaths due to ",p2,
                                  " among <5 years old for ",toupper(s2)," in ",y1), # long variable name 
                  xname="longitude",yname="latitude")
      writeRaster(ras.o5, 
                  paste0(fsav,"temp-deaths_5plus_",p1,"_",s2,"_",y1,"_",r1,"_baseline2019.nc"), #filename 
                  overwrite=T, format="CDF", datatype="FLT4S", force_v4=TRUE, compression=5,  
                  varname=paste0("mean-temp-mort_5p_",p1), # variable name
                  varunit="human deaths", 
                  longname=paste0("mean temperature-related excess deaths due to ",p2,
                                  " among 5 years & older for ",toupper(s2)," in ",y1), # long variable name 
                  xname="longitude",yname="latitude")
      
      #remove files
      rm(df1,ras.u5,ras.o5,pr.u5,pr.o5,gr1,ras.pr,temp,rval,df.u5,df.o5);gc()
      cat(paste0(i,".",j,".",k)," ")
    }
  }
}
rm(i,j,k,p1,p2,s1,s2,nm1,y1,y2,sdi1,pop1)
