######### Linear mixed models for the baseline diarrhoea projections ############

library(lme4)
library(PerformanceAnalytics)

# load data
df <- read.csv("hist_compiled_data_2000-2019_completed_w-paf_6-scalar.csv")

# create year and log outcome variables
df$yr <- df$year-2000
df$lograte <- log(df$drate+0.001) # add 0.001 for zeroes
df$age <- factor(df$age,levels=c("0-4","5-69","70+"))
p1 <- c("Shigella","Cholera","Non-typhoidal Salmonella","Enteropathogenic E coli","Enterotoxigenic E coli",
        "Campylobacter","Typhoid fever","Cryptosporidium","Rotavirus","Norovirus")
df$pathogen <- factor(df$pathogen,levels=p1)
c1 <- sort(unique(df$country))
df$country <- factor(df$country,levels=c1)
df$logscal <- log(df$scalar)
df$sdi_low <- as.numeric(df$sdi * I(df$sdi<0.8))
df$sdi_hi <- as.numeric(df$sdi * I(df$sdi>=0.8))

# build training data, aligned with CMIP6 baseline years
train.yr <- 2000:2014
train.d <- df[df$year%in%train.yr,]

# build testing data
test.yr <- 2018:2019
test.d <- df[df$year%in%test.yr,]

# base model with log scalar
mod <- lmer(formula=lograte~sdi+yr+(1|pathogen:age:country)+offset(logscal),REML=T,data=df)
coef(mod)
summary(mod)
vif.mer(mod)
kappa.mer(mod)

# check RMSE
pred1 <- predict(mod,newdata=test.d,type="response")
v1 <- exp(pred1)-0.001
sqrt(mean((test.d$drate-v1)^2))
sqrt(median((test.d$drate-v1)^2))

# save model
saveRDS(mod,"lmer_model_w-paf_6-scalar.rds")
