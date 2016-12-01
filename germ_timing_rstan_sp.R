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
  nseed<-length(unique(data$uniqueid)) #1205 unique seeds
  N<-nseed
  y<-data$daysfromstart                    # dependent variable
  temp<-data$temp   # independent variable 
  strat<-data$strat
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
  nsp<-length(unique(data$sp))
  #putting all the data together: 
  datax<-list(N=N, y=y, temp=temp, origin=origin, strat=strat,  nsp=nsp, sp=sp)
  #,nloc=nloc, nfamily=nfamily, loc=loc, family=family)
}

## fitting the stan model -------------------------------------------------

if (runstan==TRUE) {
  if (realdata==TRUE) {germdata=datax
  } else 
  {load("Fake_germdata.RData")
    germdata<-list(y=fake$y, temp=as.numeric(fake$temp), origin=as.numeric(fake$origin),
                   strat=as.numeric(fake$strat), N=nrow(fake), sp=as.numeric(fake$sp), nsp=length(unique(fake$sp)))}
  #fit_sp<- stan(file = "germdate_sp.stan", data=germdata, chains=4, iter=1200, control=list(adapt_delta=0.99))# 
  #for chains=4, iter=1200 real data has 834 divergent trans., fake has 380. for ad=.99, fake has 26, 
  #increasing the iter to 20000, with 12000 warmup, thin=2, ad=.99, get 273 divergent transisions in fake, 915 for real.  
  #fit1 <- stan(file = "germdate_sp.stan", data=germdata, chains=10, iter=1000, control = list(adapt_delta = 0.99))
  #fit_sp <- stan(file = "germdate_sp.stan", data=germdata, chains=4, iter=5000, control = list(adapt_delta = 0.99)) #high Rhat, low mixing 
  fit_sp <- stan(file = "germdate_sp.stan", data=germdata, chains=4, iter=20000, warmup=12000, thin=2,  control = list(adapt_delta = 0.99)) #This model yields 915 divergent transitions -- all below the diag in the paris plot 
  #save(fit_sp, file="germdate_sp_random3.Rdata")
  #save(fit2, file = "germdate_sp-random.Rdata")
  save(fit_sp, file="germdate_sp_fake.99_long.Rdata")
}
 
#plotting some posteriors
pdf("fit_sp_plots.pdf", width=12, height=9)
plot(fit_sp, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso","mu_b_temp", "mu_b_strat", "mu_b_origin")) 

##trace plots show poor mixing/convergence 
stan_trace(fit_sp, pars=c("mu_b_inter_to"))
stan_trace(fit_sp, pars=c("mu_b_inter_ts"))
stan_trace(fit_sp, pars=c("mu_b_inter_so"))
stan_trace(fit_sp, pars=c("mu_b_inter_tso"))
stan_trace(fit_sp, pars=c("mu_b_temp"))
stan_trace(fit_sp, pars=c("mu_b_strat"))
stan_trace(fit_sp,pars=c("mu_b_origin")) 

##pair plots show divergent transitions below the diag 

pairs(fit_sp, pars=c("mu_b_temp", "mu_b_strat", "mu_b_origin"))
pairs(fit_sp, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso"))
dev.off()
print(fit_sp, digits=3)

#----launching shiny stan---------
my_sso <- launch_shinystan(germdate_sp_random3, rstudio = getOption("shinystan.rstudio"))

