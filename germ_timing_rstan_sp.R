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
library(rstanarm)
## to download the dev version of rstanarm: 
#install.packages("devtools")
#library(devtools)
#devtools::install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)



## What do you want to do? 
runstan=TRUE      # set to true if running the stan model 
realdata=TRUE    # set to true to run on real data 

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
  temp1<-ifelse(data$temp==16.0, 1, 0) #coding temperature as binary dummy variables
  temp2<-ifelse(data$temp==20.7, 1, 0) 
  temp3<-ifelse(data$temp==25.3,1, 0) 
  strat<-ifelse(data$strat==30,0,1) 
  origin<-ifelse(data$origin=="Europe", 0, 1)
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
  loc<-as.numeric(as.factor(data$loc))
  uniqind<-as.numeric(as.factor(data$uniqind))
  nsp <- length(unique(data$sp))
  #putting all the data together: 
  datax <- data.frame(N=N, log_y=log_y, y, temp1=temp1, temp2=temp2, temp3=temp3 ,origin=origin, strat=strat,  nsp=nsp, sp=sp)
  #,nloc=nloc, nfamily=nfamily, loc=loc, family=family)
}

## fitting the stan model -------------------------------------------------

if (runstan==TRUE) {
  if (realdata==TRUE) {germdata=datax
  } else 
  {load("Fake_germdata.RData")
    germdata<-data.frame(log_y=fake$log_y, temp1=as.numeric(fake$temp1),temp2=as.numeric(fake$temp2), 
                   temp3=as.numeric(fake$temp3), origin=as.numeric(fake$origin),
                   strat=as.numeric(fake$strat), N=nrow(fake),sp=as.numeric(fake$sp),
                   nsp=length(unique(fake$sp)))}
    ##using rstanarm:
    # fitting  random intercept:
    mod_spint<-stan_lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 +
       origin*strat + origin*temp1 + origin*temp2 + origin*temp3 +
       strat*temp1 + strat*temp2 + strat*temp3 +
       origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + (1|sp), 
       data=germdata, algorithm= "sampling",
       prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5)) # by default, creates four chains with 1000 warmup, and 1000 samling 
       # note to Lizzie: Tried: prior=normal(0,30), prior_intercept=normal(0,30), prior_aux=cauchy(0,15) but no change in output
  
   #now adding random slopes
   
   mod_rs<-stan_lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 + 
                       origin*strat + origin*temp1 + origin*temp2 + origin*temp3 + 
                       strat*temp1 + strat*temp2 + strat*temp3 +
                       origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + 
                       (1|sp) +
                        (origin -1|sp) + (strat -1|sp) + (temp1 -1|sp) + (temp2 -1|sp) +  (temp3 -1|sp),
                    data=germdata, algorithm= "sampling",
                    prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),  chains=4, iter=1000)
   
   mod_rs_freq<-lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 + 
                       origin*strat + origin*temp1 + origin*temp2 + origin*temp3 + 
                       strat*temp1 + strat*temp2 + strat*temp3 +
                       origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + 
                       (1|sp) +
                       (origin -1|sp) + (strat -1|sp) + (temp1 -1|sp) + (temp2 -1|sp) +  (temp3 -1|sp),
                     data=germdata)
  
   #receive error: 
   # Error in new_CppObject_xp(fields$.module, fields$.pointer, ...) : 
   #   Exception thrown at line -1: []: accessing element out of range. index 3389 out of range; expecting index to be between 1 and 3388; index position = 1v
   # failed to create the sampler; sampling not done
   # Error in check_stanfit(stanfit) : 
   #   Invalid stanfit object produced please report bug
}


    ## see the following for running in Stan proper 
#     
#   fit_sp <- stan(file = "germdate_sp.stan", data=germdata, chains=4, iter=2000) #This model yields 915 divergent transitions -- all below the diag in the paris plot 
# 
#      #For 20000 iter, with 12000 warmup, thin=2, ad=.99, get 273 divergent transisions in fake, 915 for real.  
#   
#   #save(fit_sp, file="germdate_sp_fake.99_long.Rdata")
# }
# 
# ##pair plots show divergent transitions BELOW the diag 
# 
# pairs(fit_sp, pars=c("mu_b_temp", "mu_b_strat", "mu_b_origin"))
# pairs(fit_sp, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso"))

#----launching shiny stan---------
my_sso <- launch_shinystan(mod_rs_log, rstudio = getOption("shinystan.rstudio"))
save(mod_rs_log, file="germdate_rs_rstanarm_logged.Rdata")
mod_rs_log<-mod_rs
load("germdate_rs_rstanarm_no-log.Rdata")

#some model checking:

#jpeg(filename = "qqnorm_germdate.jpeg", height=10.5, width=8, units="in", res=500)
par(mfrow=c(2,1))
#qqnorm(resid(mod_rs_freq), main="frequentist logged")
#qqline(resid(mod_rs_freq))
qqnorm(resid(mod_rs), main= "no-log")
qqline(resid(mod_rs))
qqnorm(resid(mod_rs_log), main= "logged")
qqline(resid(mod_rs_log))
#dev.off()
