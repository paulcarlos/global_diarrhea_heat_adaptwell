# Based on Foreman et al 2018 (DOI:https://doi.org/10.1016/S0140-6736(18)31694-5)


# function 
fun.sdi <- function(values) {
  # create reference values
  rval <- data.frame("ind"=c("gdp","educ","fert"),
                     "mn"=c(250,0,3), #minimum values, GDP not in log value
                     "mx"=c(60000,17,0)) #maximum values, GDP not in log value
  
  # split values
  gdp.value <- values[1]
  educ.value <- values[2]
  fert.value <- values[3]
  
  # convert into index - GDP
  n.gdp <- log(gdp.value) - log(rval$mn[rval$ind=="gdp"]) #numerator
  d.gdp <- log(rval$mx[rval$ind=="gdp"]) - log(rval$mn[rval$ind=="gdp"])
  z.gdp <- n.gdp/d.gdp
  if (z.gdp<0.01) {
    z.gdp <- 0.01
  } else if (z.gdp>1) {
    z.gdp <- 1
  } 
  
  # convert into index - EDUC
  n.educ <- educ.value - rval$mn[rval$ind=="educ"] #numerator
  d.educ <- rval$mx[rval$ind=="educ"] - rval$mn[rval$ind=="educ"]
  z.educ <- n.educ/d.educ
  
  # convert into index - FERT
  n.fert <- fert.value - rval$mn[rval$ind=="fert"] #numerator
  d.fert <- rval$mx[rval$ind=="fert"] - rval$mn[rval$ind=="fert"]
  z.fert <- n.fert/d.fert
  
  # geometric mean
  exp(mean(log(c(z.gdp,z.educ,z.fert))))
}
