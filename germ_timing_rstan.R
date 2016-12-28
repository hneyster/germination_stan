## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())
setwd("C:/Users/Owner/Documents/GitHub/germination_stan")

##libraries
library(rstan)
library(shinystan)
library(ggplot2)

## What do you want to do? 
runstan=TRUE # set to true if running the stan model 
realdata=FALSE # set to true to run on real data 


if (realdata==TRUE) {
  
  # Setting up the  real data  for the Stan model-----------------
  
  load("germs.Rdata") #cleaned and processed real data 
  germs.y<-(subset(germs, 
                   germinated==1 & 
                     sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics 
  data<-germs.y
  N<-nrow(data)
  K<-4
  y<-data$daysfromstart                    # dependent variable
  temp<-data$temp   # independent variable 
  strat<-data$strat
  dummy_variables <- model.matrix(~ origin, data = data)
  origin<-dummy_variables[,2]
  covariates<-matrix(c(origin, temp, strat), nrow=N)                  #covariate matrix
  X = cbind(intercept=1, covariates) #covariates + intercept
  intercept<-rep(1, nrow(data))
  datax<-list(N=N, K=K, log_y=log(y), temp=temp, origin=origin, strat=strat, intercept=intercept, X=X)
}

## fitting the stan model -------------------------------------------------

if (runstan==TRUE) {
  if (realdata==TRUE) {germdata=datax
  } else 
  {load("Fake_germdata.RData")
    N=nrow(fake)
    K=4
    covariates<-matrix(c(origin=as.numeric(fake$origin), temp=as.numeric(fake$temp), strat=as.numeric(fake$strat)), nrow=N)                  #covariate matrix
    X <-  cbind(intercept=1, covariates) #covariates + intercept
    intercept<-rep(1, nrow(fake))
    germdata<-list(N=N, K=K,log_y=fake$log_y, temp=as.numeric(fake$temp), origin=as.numeric(fake$origin), 
                   strat=as.numeric(fake$strat), intercept=intercept, X=X)}
  
  fit <-stan(file = "germdate.stan", data=germdata, chains=8, iter=2000) #no divergent transitions 
}
#save(fit, file="germdate.Rdata")
## Legend:
# beta1=intercept, beta2=origin, beta3=temp, beta4=strat, beta5=temp*origin, beta6=strat*origin, beta7=strat*temp, beta8=strat*temp*origin.

#----launching shiny stan---------
my_sso <- launch_shinystan(fit, rstudio = getOption("shinystan.rstudio"))
