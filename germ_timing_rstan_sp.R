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
library(lme4)
library(ggplot2)
library(rstanarm)

source("http://peterhaschke.com/Code/multiplot.R")

## to download the dev version of rstanarm: 
#install.packages("devtools")
#library(devtools)
#devtools::install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)



## What do you want to do? 
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
  sfamily<-as.numeric(as.factor(data$uniqind))
  nsp <- length(unique(data$sp))
  #putting all the data together: 
  datax <- data.frame(N=N, log_y=log_y, y, temp1=temp1, temp2=temp2, temp3=temp3 ,origin=origin, strat=strat,  
                      nsp=nsp, sp=sp, loc=loc, sfamily=sfamily)
  #,nloc=nloc, nfamily=nfamily, loc=loc, family=family)
}

## fitting the stan model -------------------------------------------------

  if (realdata==TRUE) {
    
    germdata=datax
    
  } else 
  {load("Fake_germdata.RData")
    germdata<-data.frame(log_y=fake$log_y, temp1=as.numeric(fake$temp1),temp2=as.numeric(fake$temp2), 
                   temp3=as.numeric(fake$temp3), origin=as.numeric(fake$origin),
                   strat=as.numeric(fake$strat), N=nrow(fake),sp=as.numeric(fake$sp),
                   nsp=length(unique(fake$sp)))}
    ## Fitting models with rstanarm:
  
  
    # Model 1: random intercept:
    mod__time_spint<-stan_lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 +
       origin*strat + origin*temp1 + origin*temp2 + origin*temp3 +
       strat*temp1 + strat*temp2 + strat*temp3 +
       origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + (1|sp), 
       data=germdata, algorithm= "sampling",
       prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5)) # by default, creates four chains with 1000 warmup, and 1000 samling 
       # note to Lizzie: Tried: prior=normal(0,30), prior_intercept=normal(0,30), prior_aux=cauchy(0,15) but no change in output
  
   #Model 2: now adding random slopes
    
    mod_time_rslope<-stan_lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 + 
                        origin*strat + origin*temp1 + origin*temp2 + origin*temp3 + 
                        strat*temp1 + strat*temp2 + strat*temp3 +
                        origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + 
                        (1|sp) +
                        (origin -1|sp) + (strat -1|sp) + (temp1 -1|sp) + (temp2 -1|sp) +  (temp3 -1|sp) + (origin:strat -1|sp) + (origin:temp1 -1|sp) + (origin:temp2 -1|sp) + (origin:temp3 -1|sp)+
                        (strat:temp1 -1|sp) + (strat:temp2 -1|sp) + (strat:temp3 -1|sp) + (origin:strat:temp1 -1|sp) + (origin:strat:temp2 -1|sp) + (origin:strat:temp3 -1|sp),
                      data=germdata, algorithm= "sampling",
                      prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),  chains=4, iter=2000)
   
  # Model 3: now adding additional nested random effects: 
    
   mod_time_log<-stan_lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 + 
                            origin:strat + origin:temp1 + origin:temp2 + origin:temp3 + 
                            strat:temp1 + strat:temp2 + strat:temp3 +
                            origin:strat:temp1 +  origin:strat:temp2 + origin:strat:temp3 + 
                            (1|sp/loc/sfamily) +
                            (origin -1|sp/loc/sfamily) + (strat -1|sp/loc/sfamily) + (temp1 -1|sp/loc/sfamily) + 
                            (temp2 -1|sp/loc/sfamily) +  (temp3 -1|sp/loc/sfamily)+
                              (origin:strat -1|sp/loc/sfamily) + (origin:temp1 -1|sp/loc/sfamily) + (origin:temp2 -1|sp/loc/sfamily) + 
                              (origin:temp3 -1|sp/loc/sfamily) +  (strat:temp1 -1|sp/loc/sfamily) + (strat:temp2 -1|sp/loc/sfamily) + 
                              (strat:temp3 -1|sp/loc/sfamily) + (origin:strat:temp1 -1|sp/loc/sfamily) + 
                              (origin:strat:temp2 -1|sp/loc/sfamily) + (origin:strat:temp3 -1|sp/loc/sfamily),
                          data=germdata, algorithm= "sampling",
                    prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),  chains=4, iter=2000)
   
   # model 4: now trying Poisson error distribution: 
   
   mod_time_pois<-stan_glmer(y ~ origin + strat + temp1 + temp2 + temp3 + 
                               origin:strat + origin:temp1 + origin:temp2 + origin:temp3 + 
                               strat:temp1 + strat:temp2 + strat:temp3 +
                               origin:strat:temp1 +  origin:strat:temp2 + origin:strat:temp3 + 
                               (1|sp/loc/sfamily) +
                               (origin -1|sp/loc/sfamily) + (strat -1|sp/loc/sfamily) + (temp1 -1|sp/loc/sfamily) + 
                               (temp2 -1|sp/loc/sfamily) +  (temp3 -1|sp/loc/sfamily)+
                               (origin:strat -1|sp/loc/sfamily) + (origin:temp1 -1|sp/loc/sfamily) + (origin:temp2 -1|sp/loc/sfamily) + 
                               (origin:temp3 -1|sp/loc/sfamily) +  (strat:temp1 -1|sp/loc/sfamily) + (strat:temp2 -1|sp/loc/sfamily) + 
                               (strat:temp3 -1|sp/loc/sfamily) + (origin:strat:temp1 -1|sp/loc/sfamily) + 
                               (origin:strat:temp2 -1|sp/loc/sfamily) + (origin:strat:temp3 -1|sp/loc/sfamily),
                             data=germdata, algorithm= "sampling", family=poisson,
                             prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),  chains=4, iter=2000)
   
  #and checking against lmer model:
   mod_rs_freq2<-lmer(log_y ~ origin + strat + temp1 + temp2 + temp3 + 
                                  origin*strat + origin*temp1 + origin*temp2 + origin*temp3 + 
                                  strat*temp1 + strat*temp2 + strat*temp3 +
                                  origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + 
                                  (1|sp) +
                                  (origin -1|sp) + (strat -1|sp) + (temp1 -1|sp) + (temp2 -1|sp) +  (temp3 -1|sp) + (origin:strat -1|sp) + (origin:temp1 -1|sp) + (origin:temp2 -1|sp) + (origin:temp3 -1|sp)+
                                  (strat:temp1 -1|sp) + (strat:temp2 -1|sp) + (strat:temp3 -1|sp) + (origin:strat:temp1 -1|sp) + (origin:strat:temp2 -1|sp) + (origin:strat:temp3 -1|sp),
                                data=germdata, algorithm= "sampling",
                                prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),  chains=4, iter=2000)

   save(mod_time_log, "mod_time_log.Rdata")

#------### MODEL ANALYSIS ###-------------
   
#launching shiny stan
   
load("mod_time_log.Rdata")
my_sso <- launch_shinystan(mod_rs, rstudio = getOption("shinystan.rstudio"))

mod_to_plot<- mod_time_pois

p1<-plot(mod_to_plot, pars=c("origin", "strat", "temp1", "temp2", "temp3", "origin:strat", "origin:temp1", "origin:temp2",
                             "origin:temp3", "strat:temp1", "strat:temp2", "strat:temp3", "origin:strat:temp1", "origin:strat:temp2", "origin:strat:temp3")
p2<-plot(mod_to_plot, pars=c("b[origin sp:1]", "b[strat sp:1]", "b[temp1 sp:1]", "b[temp2 sp:1]", "b[temp3 sp:1]", "b[origin:strat sp:1]",
                             "b[origin:temp1 sp:1]", "b[origin:temp2 sp:1]", "b[origin:temp3 sp:1]", "b[strat:temp1 sp:1]", "b[strat:temp2 sp:1]", "b[strat:temp3 sp:1]", "b[origin:strat:temp1 sp:1]",
                             "b[origin:strat:temp2 sp:1]", "b[origin:strat:temp3 sp:1]"))
p3<-plot(mod_to_plot, pars=c("b[origin sp:2]", "b[strat sp:2]", "b[temp1 sp:2]", "b[temp2 sp:2]", "b[temp3 sp:2]", "b[origin:strat sp:2]",
                             "b[origin:temp1 sp:2]", "b[origin:temp2 sp:2]", "b[origin:temp3 sp:2]","b[strat:temp1 sp:2]", "b[strat:temp2 sp:2]", "b[strat:temp3 sp:2]", "b[origin:strat:temp1 sp:2]", 
                             "b[origin:strat:temp2 sp:2]", "b[origin:strat:temp3 sp:2]"))
p4<-plot(mod_to_plot, pars=c("b[origin sp:3]", "b[strat sp:3]", "b[temp1 sp:3]", "b[temp2 sp:3]", "b[temp3 sp:3]", "b[origin:strat sp:3]",
                             "b[origin:temp1 sp:3]", "b[origin:temp2 sp:3]",  "b[origin:temp3 sp:3]", "b[strat:temp1 sp:3]", "b[strat:temp2 sp:3]", "b[strat:temp3 sp:3]", "b[origin:strat:temp1 sp:3]",
                             "b[origin:strat:temp2 sp:3]", "b[origin:strat:temp3 sp:3]"))
p5<-plot(mod_to_plot, pars=c("b[origin sp:4]", "b[strat sp:4]", "b[temp1 sp:4]", "b[temp2 sp:4]", "b[temp3 sp:4]", "b[origin:strat sp:4]",
                             "b[origin:temp1 sp:4]", "b[origin:temp2 sp:4]",  "b[origin:temp3 sp:4]", "b[strat:temp1 sp:4]", "b[strat:temp2 sp:4]", "b[strat:temp3 sp:4]", "b[origin:strat:temp1 sp:4]",
                             "b[origin:strat:temp2 sp:4]", "b[origin:strat:temp3 sp:4]"))
p6<-plot(mod_to_plot, pars=c("b[origin sp:5]", "b[strat sp:5]", "b[temp1 sp:5]", "b[temp2 sp:5]", "b[temp3 sp:5]", "b[origin:strat sp:5]",
                             "b[origin:temp1 sp:5]", "b[origin:temp2 sp:5]",  "b[origin:temp3 sp:5]", "b[strat:temp1 sp:5]", "b[strat:temp2 sp:5]", "b[strat:temp3 sp:5]", "b[origin:strat:temp1 sp:5]",
                             "b[origin:strat:temp2 sp:5]", "b[origin:strat:temp3 sp:5]"))
p7<-plot(mod_to_plot, pars=c("b[origin sp:6]", "b[strat sp:6]", "b[temp1 sp:6]", "b[temp2 sp:6]", "b[temp3 sp:6]", "b[origin:strat sp:6]",
                             "b[origin:temp1 sp:6]", "b[origin:temp2 sp:6]",  "b[origin:temp3 sp:6]","b[strat:temp1 sp:6]", "b[strat:temp2 sp:6]", "b[strat:temp3 sp:6]", "b[origin:strat:temp1 sp:6]",
                             "b[origin:strat:temp2 sp:6]", "b[origin:strat:temp3 sp:6]"))
p8<-plot(mod_to_plot, pars=c("b[origin sp:7]", "b[strat sp:7]", "b[temp1 sp:7]", "b[temp2 sp:7]", "b[temp3 sp:7]", "b[origin:strat sp:7]",
                             "b[origin:temp1 sp:7]", "b[origin:temp2 sp:7]",  "b[origin:temp3 sp:7]","b[strat:temp1 sp:7]", "b[strat:temp2 sp:7]", "b[strat:temp3 sp:7]", "b[origin:strat:temp1 sp:7]",
                             "b[origin:strat:temp2 sp:7]", "b[origin:strat:temp3 sp:7]"))

pdf("timing_logged_rsi.pdf", width=20, height=11)
multiplot(p1,p2, p3, p4, p5, p6, p7, p8 cols=4)
dev.off()


#some model checking:
posterior_interval(mod_time_pois, prob = 0.95, type = "central")
#jpeg(filename = "qqnorm_germdate.jpeg", height=10.5, width=8, units="in", res=500)
par(mfrow=c(2,1))
#qqnorm(resid(mod_rs_freq), main="frequentist logged")
#qqline(resid(mod_rs_freq))
qqnorm(resid(mod_time), main= "no-log")
qqline(resid(mod_time))
qqnorm(resid(mod_time_log), main= "logged")
qqline(resid(mod_time_log))
#dev.off()
