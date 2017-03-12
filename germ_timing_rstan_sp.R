## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/undergrads/harold/analyses/germination_stan") 
} else 
  setwd("C:/Users/Owner/Documents/GitHub/germination_stan")

##libraries
library(rstan)
library(shinystan)
library(ggplot2)

## What do you want to do? 
runstan=FALSE      # set to true if running the stan model 
realdata=FALSE    # set to true to run on real data 

if (realdata==TRUE) {
  # Setting up the  real data  for the Stan model-----------------
  
  load("germs.Rdata") #cleaned and processed real data 
  germs.y<-(subset(germs, 
                   germinated==1 & 
                     sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics 
  data<-germs.y
  nseed<-length(unique(data$uniqueid)) #1205 unique seeds
  N<-nseed
  y<-data$daysfromstart    # dependent variable
  log_y=log(y)
  temp<-data$temp   
  dx$chill1 = ifelse(dx$chill == 2, 1, 0)
  dx$chill2 = ifelse(dx$chill == 3, 1, 0)
  strat= data$strat
  strat= ifelse(strat==30,1,2)
  dummy_variables <- model.matrix(~ origin, data = data) 
  origin<-dummy_variables[,2]
  intercept<-rep(1, nrow(data))
  #setting up to random effects data:
  nsp<-length(unique(germs.y$sp))
  sp_alph<-data$sp
  
  sp<-ifelse (sp_alph=="CAPBUR", 1,     #making sp numeric, in alphabetical order 
              ifelse(sp_alph=="CHEMAJ",2,
                     ifelse(sp_alph=="DACGLO", 3, 
                            ifelse(sp_alph=="PLALAN", 4,
                                   ifelse(sp_alph=="PLAMAJ", 5, 
                                          ifelse(sp_alph=="RUMCRI", 6, 7))))))
  nsp <- length(unique(data$sp))
  #putting all the data together: 
  datax <- list(N=N, log_y=log_y, temp=temp, origin=origin, strat=strat,  nsp=nsp, sp=sp)
  #,nloc=nloc, nfamily=nfamily, loc=loc, family=family)
}

## fitting the stan model -------------------------------------------------

if (runstan==TRUE) {
  if (realdata==TRUE) {germdata=datax
  } else 
  {load("Fake_germdata.RData")
    germdata<-list(log_y=fake$log_y, temp=as.numeric(fake$temp), origin=as.numeric(fake$origin),
                   strat=as.numeric(fake$strat), N=nrow(fake), sp=as.numeric(fake$sp), nsp=length(unique(fake$sp)))}
  
  fit_sp <- stan(file = "germdate_sp.stan", data=germdata, chains=4, iter=2000) #This model yields 915 divergent transitions -- all below the diag in the paris plot 

     #For 20000 iter, with 12000 warmup, thin=2, ad=.99, get 273 divergent transisions in fake, 915 for real.  
  
  #save(fit_sp, file="germdate_sp_fake.99_long.Rdata")
}

##pair plots show divergent transitions BELOW the diag 

pairs(fit_sp, pars=c("mu_b_temp", "mu_b_strat", "mu_b_origin"))
pairs(fit_sp, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso"))


#----launching shiny stan---------
my_sso <- launch_shinystan(fit_sp, rstudio = getOption("shinystan.rstudio"))

