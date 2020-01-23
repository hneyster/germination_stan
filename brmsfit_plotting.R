#Code to plot model coefficients for germination time with a brms.fit object. 
#Written by Harold Eyster 1/15/18

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())

setwd("C:/Users/Owner/Documents/GitHub/germination_stan")
source("http://peterhaschke.com/Code/multiplot.R") #so that the multiplot function works 
source("https://raw.githubusercontent.com/jaredlander/coefplot/master/R/position.r") # for vertical dodging 



##libraries
library(rstan)
library(lme4)
library(ggplot2)
library(rstanarm)
library(brms)
library(tidyr)
library(RCurl)
library(forcats)



####### Days from Germination Plot:
#######
#######

#data
load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois_brm.Rdata") #this is my brms stanfit object. 

m<-mod_time_pois_brm 
sum.m<-summary(m) #summary of the object 
cri.f<-as.data.frame(sum.m$fixed[,c("Estimate", "l-95% CI", "u-95% CI")]) #Extracting the global effects with CI
cri.f<-cri.f[-1,] #removing the intercept 
fdf1<-as.data.frame(rbind(as.vector(cri.f[,1]), as.vector(cri.f[,2]), as.vector(cri.f[,3]))) # Making into a dataframe 
fdf2<-cbind(fdf1, c(0, 0, 0) , c("Estimate", "2.5%", "97.5%")) #adding zeros to denote that these are the global effects 
names(fdf2)<-c(rownames(cri.f), "sp", "perc") #renaming 

cri.r<-(ranef(m, summary = TRUE, robust = FALSE,
              probs = c(0.025, 0.975)))$sp   #extracting the random effects 
cri.r2<-cri.r[, ,-1] #removing the intercept 
cri.r2<-cri.r2[,-2,] #removing the est error
dims<-dim(cri.r2)
twoDimMat <- matrix(cri.r2, prod(dims[1:2]), dims[3]) #converting to matrix 
mat2<-cbind(twoDimMat, c(rep(1:7, length.out=21)), rep(c("Estimate", "2.5%", "97.5%"), each=7)) #adding a species column 
df<-as.data.frame(mat2)
names(df)<-c(rownames(cri.f), "sp", "perc") #renaming 
dftot<-rbind(fdf2, df) #combining  fixed effects and random effects 
dflong<- gather(dftot, var, value, origin:`origin:strat:temp3`, factor_key=TRUE) #converting to long format

#adding the fixed effect  estiamtes to the random effect values 
for (i in seq(from=1,to=nrow(dflong), by=24)) {
  for (j in seq(from=3, to=23, by=1)) {
    dflong$value[i+j]<- as.numeric(dflong$value[i+j]) + as.numeric(dflong$value[i])
  }
}
dflong$rndm<-ifelse(dflong$sp>0, 2, 1) #adding a new column that signifies weather a coeff is a fixed or a random effect 
dfwide<-spread(dflong, perc, value) #widening it a bit
dfwide[,4:6] <- as.data.frame(lapply(c(dfwide[,4:6]), as.numeric )) #making numeric 
dfwide$sp<-as.factor(dfwide$sp) 
## plotting

pd <- position_dodgev(height = -0.5) #the negative here makes it so the fixed effects are plotted above the others 


fig1 <-ggplot(dfwide, aes(x=Estimate, y=var, color=factor(sp), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels = c("Fixed effects", "CAPBUR", "CHEMAJ", "DACGLO", "PLALAN", "PLAMAJ", "RUMCRI", "TAROFF"),
                      values=c("blue", "red", "orangered1", "sienna4", "green4", "green1", "purple2", "magenta2"))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Days to Germination")+ 
  scale_y_discrete(limits = rev(unique(sort(dfwide$var)))) #reverses the y axis 
#pdf(file.path(figpath, "Fig1.pdf"), width = 7, height = 8)
fig1
#dev.off()

#####
#####
#####Germination Rate Plot 
