## code to run the APC function and produce tables ## 

rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(arm)
library(rstan)
library(rstanarm)
library(reshape2)
library(ggplot2)
library(xtable)

## load function ## 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("germ_apc.R")

##################################
#### APCs for Growth Rate #### 
#################################
apc_subset<-function(sub){ 
    load("C:/Users/Owner/Documents/Thesis/Stan/mod_gr.Rdata")
    load("C:/Users/Owner/Documents/github/germination_stan/datax.Rdata") 
    datax<-datax[datax$origin==sub,]
    gr_s<-apc(mod_gr,"strat",datax,type="binary")
    gr_t1<-apc(mod_gr,"temp1", datax, type="binary")
    gr_t2<-apc(mod_gr,"temp2", datax, type="binary")
    gr_t3<-apc(mod_gr,"temp3", datax, type="binary")
    gr_sp<-apc(mod_gr, "sp",datax,type="categorical")
  
  ###################################
  ### APCs for germ rate model######## 
  ###################################
  
  
  load("C:/Users/Owner/Documents/Thesis/Stan/mod_rate.Rdata") 
  load("C:/Users/Owner/Documents/github/germination_stan/rate_data.Rdata") 
  rate_data<-rate_data[rate_data$origin==sub,]
  ge_s<-apc(mod_rate,"strat",rate_data,type="binary")
  ge_sp<-apc(mod_rate,"sp",rate_data,type="categorical")
  ge_t1<-apc(mod_rate,"temp1", rate_data,type="binary")
  ge_t2<-apc(mod_rate,"temp2", rate_data,type="binary")
  ge_t3<-apc(mod_rate,"temp3", rate_data,type="binary")
  
  ###################################
  #### APCs for germ timing model ###
  ###################################
  
  load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois.Rdata") #this is a brms stanfit object. 
  load("C:/Users/Owner/Documents/github/germination_stan/time_data.rdata") 
  time_data<-time_data[time_data$origin==sub,]
  t_s<-apc(mod_time_pois,"strat",time_data,type = "binary")
  t_sp<-apc(mod_time_pois, "sp", time_data, type="categorical")
  t_t1<-apc(mod_time_pois,"temp1", time_data, type = "binary")
  t_t2<-apc(mod_time_pois,"temp2", time_data, type="binary")
  t_t3<-apc(mod_time_pois,"temp3", time_data, type="binary")
  
  
  #################################
  ### Making a table #############
  ################################
  apctable<-rbind(cbind(ge_s,ge_t1, ge_t2, ge_t3,ge_sp),
                  cbind(t_s,t_t1,t_t2,t_t3,t_sp),
                  cbind( gr_s, gr_t1,gr_t2,gr_t3,gr_sp))
  apctable<-round(apctable,3)
  table1<-as.data.frame(rbind(
    cbind(paste(apctable[1,1],"$\\pm$",apctable[1,2]),
          paste(apctable[1,3], "$\\pm$",apctable[1,4]),
          paste(apctable[1,5], "$\\pm$",apctable[1,6]),
          paste(apctable[1,7], "$\\pm$",apctable[1,8]),                                            
          paste(apctable[1,9], "$\\pm$",apctable[1,10])),
    cbind(paste(apctable[2,1],"$\\pm$",apctable[2,2]),
          paste(apctable[2,3], "$\\pm$",apctable[2,4]),
          paste(apctable[2,5], "$\\pm$",apctable[2,6]),
          paste(apctable[2,7], "$\\pm$",apctable[2,8]),                                            
          paste(apctable[2,9], "$\\pm$",apctable[2,10])),
    cbind(paste(apctable[3,1],"$\\pm$",apctable[3,2]),
          paste(apctable[3,3], "$\\pm$",apctable[3,4]),
          paste(apctable[3,5], "$\\pm$",apctable[3,6]),
          paste(apctable[3,7], "$\\pm$",apctable[3,8]),                                            
          paste(apctable[3,9], "$\\pm$",apctable[3,10]))
  ))
  names(table1)<-c( "stratification", "temperature 1", "temperature 2","temperature 3","species")
  row.names(table1)<-c("germination rate (\%)", "germination date (days)", "growth rate (cm/day)")
  table_t<-t(table1)
  return(print(xtable(table_t),tabular.environment = "longtable",include.rownames = TRUE, floating = FALSE, 
        sanitize.text.function = identity))
}

apc_subset(sub=0)
apc_subset(sub=1)
setwd("C:/Users/Owner/Documents/github/germination_stan")


#############################
######### SAVING ###########
############################

list<-ls()
list<-list[-which(list %in% c("mod_gr","mod_rate","mod_time_pois"))] # list all variables except the one you want to save
save(list = list, file = "apc_germ.RData")

########################################################
#### NOW SUBSETTING ORIGIN TO TRY AND GET A LOOK AT PLASTICITY #### 
###################################################

## first europe ## 

