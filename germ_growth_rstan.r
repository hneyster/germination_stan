## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())

##libraries
library(rstan)
library(shinystan)
library(lme4)
library(ggplot2)
library(rstanarm)
library(grid)
library(here)
#source("http://peterhaschke.com/Code/multiplot.R")
## to download the dev version of rstanarm: 
#install.packages("devtools")
#library(devtools)
#devtools::install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)



## What do you want to do? 
runstan=TRUE      # set to true if running the stan model 
realdata=TRUE    # set to true to run on real data 

if (realdata==TRUE) {
  # Setting up the  real data  for the Stan model-----------------
  
  load("hlms.Rdata")  #cleaned and processed real height data 
  hlms$uniqind<- do.call(paste, c(hlms[c("sp", "individual")], sep="_")) #creates a column that defines the individual 
  
  data<-hlms
  nseed<-length(unique(data$uniqueid)) #1578 unique seeds
  N<-nseed
  y<-data$gr  # dependent variable
  temp1<-ifelse(data$temp==16.0, 1, 0) #coding temperature as binary dummy variables
  temp2<-ifelse(data$temp==20.7, 1, 0) 
  temp3<-ifelse(data$temp==25.3,1, 0) 
  strat<-ifelse(data$strat==60,0,1) # long strat/winter is the reference level (changed 2021-04-18)
  origin<-ifelse(data$origin=="Europe", 0, 1)
  intercept<-rep(1, nrow(data))
  #setting up to random effects data:
  nsp<-length(unique(data$sp))
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
  datax <- data.frame(N=N, y=y, temp1=temp1, temp2=temp2, temp3=temp3 ,origin=origin, strat=strat,  
                      nsp=nsp, sp=sp, loc=loc, sfamily=sfamily)
 # save(datax,file="datax.rdata")
}

## fitting the stan model -------------------------------------------------

if (runstan==TRUE) {
  germdata=datax
  
  ##using rstanarm:
  #first using no just random intercept
  
  # fitting  random intercept:
  mod_gr.i<-stan_lmer(y ~ origin + strat + temp1 + temp2 + temp3 +
                          origin*strat + origin*temp1 + origin*temp2 + origin*temp3 +
                          strat*temp1 + strat*temp2 + strat*temp3 +
                          origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 + (1|sp), 
                        data=germdata, algorithm= "sampling",
                        prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),
                        chains=1, iter=200) # by default, creates four chains with 1000 warmup, and 1000 samling 
  
  
  #now adding random slopes
  # this is model used in paper: 
  
  mod_gr<-stan_lmer(y ~ origin + strat + temp1 + temp2 + temp3 + 
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
                        prior=normal(), prior_intercept=normal(0,10), prior_aux=cauchy(0,5),
                        chains=4, iter=2000) # by default, creates four chains with 1000 warmup, and 1000 samling 
            #caused one divergent transition after warm-up
  save(mod_gr, file="mod_gr.rdata")
  
          
  #and checking against lmer model:
  mod_grfreq2<-lmer(y ~ origin + strat + temp1 + temp2 + temp3 + 
                       origin:strat + origin:temp1 + origin:temp2 + origin:temp3 + 
                       strat:temp1 + strat:temp2 + strat:temp3 +
                       origin:strat:temp1 +  origin:strat:temp2 + origin:strat:temp3 + 
                       (1|sp/loc/sfamily) +
                       (origin -1|sp/loc/sfamily) + (strat -1|sp/loc/sfamily) + (temp1 -1|sp/loc/sfamily) + 
                       (temp2 -1|sp/loc/sfamily) +  (temp3 -1|sp/loc/sfamily),
                     data=germdata)
  
}

  #----launching shiny stan---------
  load("mod_gr.Rdata")
  my_sso <- launch_shinystan(mod_gr, rstudio = getOption("shinystan.rstudio"))
  
  
  #plotting
  p1<-plot(mod_gr, pars=c("origin", "strat", "temp1", "temp2", "temp3", "origin:strat", "origin:temp1", "origin:temp2",
                            "origin:temp3", "strat:temp1", "strat:temp2", "strat:temp3", "origin:strat:temp1", "origin:strat:temp2", "origin:strat:temp3"))
  
  p2<-plot(mod_gr, pars=c("b[origin sp:1]", "b[strat sp:1]", "b[temp1 sp:1]", "b[temp2 sp:1]", "b[temp3 sp:1]", "b[origin:strat sp:1]",
                            "b[origin:temp1 sp:1]", "b[origin:temp2 sp:1]", "b[origin:temp3 sp:1]", "b[strat:temp1 sp:1]", "b[strat:temp2 sp:1]", "b[strat:temp3 sp:1]", "b[origin:strat:temp1 sp:1]",
                            "b[origin:strat:temp2 sp:1]", "b[origin:strat:temp3 sp:1]"))
  
  p3<-plot(mod_gr, pars=c("b[origin sp:2]", "b[strat sp:2]", "b[temp1 sp:2]", "b[temp2 sp:2]", "b[temp3 sp:2]", "b[origin:strat sp:2]",
                            "b[origin:temp1 sp:2]", "b[origin:temp2 sp:2]", "b[origin:temp3 sp:2]","b[strat:temp1 sp:2]", "b[strat:temp2 sp:2]", "b[strat:temp3 sp:2]", "b[origin:strat:temp1 sp:2]", 
                            "b[origin:strat:temp2 sp:2]", "b[origin:strat:temp3 sp:2]"))
  
  p4<-plot(mod_gr, pars=c("b[origin sp:3]", "b[strat sp:3]", "b[temp1 sp:3]", "b[temp2 sp:3]", "b[temp3 sp:3]", "b[origin:strat sp:3]",
                            "b[origin:temp1 sp:3]", "b[origin:temp2 sp:3]",  "b[origin:temp3 sp:3]", "b[strat:temp1 sp:3]", "b[strat:temp2 sp:3]", "b[strat:temp3 sp:3]", "b[origin:strat:temp1 sp:3]",
                            "b[origin:strat:temp2 sp:3]", "b[origin:strat:temp3 sp:3]"))
  
  p5<-plot(mod_gr, pars=c("b[origin sp:4]", "b[strat sp:4]", "b[temp1 sp:4]", "b[temp2 sp:4]", "b[temp3 sp:4]", "b[origin:strat sp:4]",
                            "b[origin:temp1 sp:4]", "b[origin:temp2 sp:4]",  "b[origin:temp3 sp:4]", "b[strat:temp1 sp:4]", "b[strat:temp2 sp:4]", "b[strat:temp3 sp:4]", "b[origin:strat:temp1 sp:4]",
                            "b[origin:strat:temp2 sp:4]", "b[origin:strat:temp3 sp:4]"))
  
  p6<-plot(mod_gr, pars=c("b[origin sp:5]", "b[strat sp:5]", "b[temp1 sp:5]", "b[temp2 sp:5]", "b[temp3 sp:5]", "b[origin:strat sp:5]",
                            "b[origin:temp1 sp:5]", "b[origin:temp2 sp:5]",  "b[origin:temp3 sp:5]", "b[strat:temp1 sp:5]", "b[strat:temp2 sp:5]", "b[strat:temp3 sp:5]", "b[origin:strat:temp1 sp:5]",
                            "b[origin:strat:temp2 sp:5]", "b[origin:strat:temp3 sp:5]"))
  
  p7<-plot(mod_gr, pars=c("b[origin sp:6]", "b[strat sp:6]", "b[temp1 sp:6]", "b[temp2 sp:6]", "b[temp3 sp:6]", "b[origin:strat sp:6]",
                            "b[origin:temp1 sp:6]", "b[origin:temp2 sp:6]",  "b[origin:temp3 sp:6]","b[strat:temp1 sp:6]", "b[strat:temp2 sp:6]", "b[strat:temp3 sp:6]", "b[origin:strat:temp1 sp:6]",
                            "b[origin:strat:temp2 sp:6]", "b[origin:strat:temp3 sp:6]"))
  # pdf("gr.pdf", width=20, height=11)
  multiplot(p1,p2, p3, p4, p5, p6, p7, cols=4)
  #dev.off()
  #some model checking:
  plot(mod_rs_log)
  plot(mod_rs)
  #jpeg(filename = "qqnorm_germdate.jpeg", height=10.5, width=8, units="in", res=500)
  par(mfrow=c(2,1))
  #qqnorm(resid(mod_rs_freq), main="frequentist logged")
  #qqline(resid(mod_rs_freq))
  qqnorm(resid(mod_rs), main= "no-log")
  qqline(resid(mod_rs))
  qqnorm(resid(mod_rs_log), main= "logged")
  qqline(resid(mod_rs_log))
  #dev.off()
  