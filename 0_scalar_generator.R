############### relative risks for risk factors #################
# create dataframe
rr.hi <- data.frame("wat.unimp"=11.084,"sani.unimp"=3.203,"no.hw"=1.522,"stunt"=1.222,"ors"=1/0.31,
                    "rotavac"=1/0.0569)
rr.mi <- data.frame("wat.unimp"=11.084,"sani.unimp"=3.203,"no.hw"=1.522,"stunt"=1.222,"ors"=1/0.31,
                    "rotavac"=1/0.2398)
rr.li <- data.frame("wat.unimp"=11.084,"sani.unimp"=3.203,"no.hw"=1.522,"stunt"=1.222,"ors"=1/0.31,
                    "rotavac"=1/0.5525)
# RR for ORS based on Munos 2010 using 69% relative reduction
# RR for rotavirus vaccination is from Sun et al 2021, using RCTs 1 year follow-up
# RR for WASH and stunting from Murray et al 2020


# income categories 
wbinc <- read.csv("countries_wb-category_1987-2022.csv")

# function for rota
yr=2019;ctry="PHL";wat=0;san=0.2;hw=0.1;stu=0.4;or=0.04;rvac=0.2
fu5rota <- function(yr,ctry,wat,san,hw,stu,or,rvac) {
  inc <- wbinc$wbcat[wbinc$iso3==ctry & wbinc$year==yr & !is.na(wbinc$wbcat)]
  if (inc=="H") {
    rval <- rr.hi
  } else if (inc=="L") {
    rval <- rr.li
  } else {
    rval <- rr.mi
  }
  wat=as.numeric(wat)*(rval$wat.unimp-1); san=as.numeric(san)*(rval$sani.unimp-1); 
  hw=as.numeric(hw)*(rval$no.hw-1); stu=as.numeric(stu)*(rval$stunt-1); or=as.numeric(or)*(rval$ors-1); 
  rvac=as.numeric(rvac)*(rval$rotavac-1)
  vwat <- wat/(wat+1)
  vsan <- san/(san+1)
  vhw <- hw/(hw+1)
  vstu <- stu/(stu+1)
  vor <- or/(or+1)
  vrvac <- rvac/(rvac+1)
  v1 <- c(vwat,vsan,vhw,vstu,vor,vrvac)
  if (all(is.na(v1))) {
    scal <- 0
  } else {
    paf1 <- 1-(prod(1-v1[which(v1>0)]))
    scal <- log(1/(1-paf1))
  }
  scal
}

fu5all <- function(wat,san,hw,stu,or) {
  rval <- rr.mi
  wat=as.numeric(wat)*(rval$wat.unimp-1); san=as.numeric(san)*(rval$sani.unimp-1); 
  hw=as.numeric(hw)*(rval$no.hw-1); stu=as.numeric(stu)*(rval$stunt-1); or=as.numeric(or)*(rval$ors-1); 
  vwat <- wat/(wat+1)
  vsan <- san/(san+1)
  vhw <- hw/(hw+1)
  vstu <- stu/(stu+1)
  vor <- or/(or+1)
  v1 <- c(vwat,vsan,vhw,vstu,vor)
  if (all(is.na(v1))) {
    scal <- 0
  } else {
    paf1 <- 1-(prod(1-v1[which(v1>0)]))
    scal <- log(1/(1-paf1))
  }
  scal
}

fo5 <- function(wat,san,hw) {
  rval <- rr.mi
  wat=as.numeric(wat)*(rval$wat.unimp-1); san=as.numeric(san)*(rval$sani.unimp-1); 
  hw=as.numeric(hw)*(rval$no.hw-1) 
  vwat <- wat/(wat+1)
  vsan <- san/(san+1)
  vhw <- hw/(hw+1)
  v1 <- c(vwat,vsan,vhw)
  if (all(is.na(v1))) {
    scal <- 0
  } else {
    paf1 <- 1-(prod(1-v1[which(v1>0)]))
    scal <- log(1/(1-paf1))
  }
  scal
}
