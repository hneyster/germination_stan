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
source("apc_withinsp.R")
##################################
#### APCs for Growth Rate #### 
#################################
load("C:/Users/Owner/Documents/Thesis/Stan/mod_gr.Rdata")
load("C:/Users/Owner/Documents/github/germination_stan/datax.Rdata") 
gr_o<-apc(mod_gr,"origin",datax,type = "binary")
gr_s<-apc(mod_gr,"strat",datax,type="binary")
gr_sp<-apc(mod_gr, "sp",datax,type="categorical", nested=c("loc", "sfamily"))
gr_pop <-apc(mod_gr, "loc", datax, type = "categorical", nested = c("sp", "sfamily"))
gr_t1<-apc(mod_gr,"temp1", datax, type="binary")
gr_t2<-apc(mod_gr,"temp2", datax, type="binary")
gr_t3<-apc(mod_gr,"temp3", datax, type="binary")

## Just for PLALAN:
gr_plalan <- apc_withinsp(mod_gr, "loc", datax,type = "categorical", nested="sfamily", sp=4)

## now looping across all species: 
gr_withinsp_apc<-data.frame()
for (i in 1:7){
a <- apc_withinsp(mod_gr, "loc", datax,type = "categorical", nested="sfamily", sp=i)
gr_withinsp_apc = rbind(gr_withinsp_apc, a)
}

###################################
### APCs for germ rate model######## 
###################################


load("C:/Users/Owner/Documents/Thesis/Stan/mod_rate.Rdata") 
load("C:/Users/Owner/Documents/github/germination_stan/rate_data.Rdata") 
ge_o<-apc(mod_rate,"origin",rate_data,type="binary")
ge_s<-apc(mod_rate,"strat",rate_data,type="binary")
ge_sp<-apc(mod_rate,"sp",rate_data,type="categorical", nested=c("loc", "sfamily"))
ge_pop <-apc(mod_rate, "loc", rate_data, type = "categorical", nested = c("sp", "sfamily"))
ge_t1<-apc(mod_rate,"temp1", rate_data,type="binary")
ge_t2<-apc(mod_rate,"temp2", rate_data,type="binary")
ge_t3<-apc(mod_rate,"temp3", rate_data,type="binary")

ge_plalan<-apc_withinsp(mod_rate, "loc", rate_data,type = "categorical", nested="sfamily", sp=4)

###################################
#### APCs for germ timing model ###
###################################

load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois.Rdata") 
load("C:/Users/Owner/Documents/github/germination_stan/time_data.rdata") 
t_o<-apc(mod_time_pois,"origin",time_data,type = "binary")
t_s<-apc(mod_time_pois,"strat",time_data,type = "binary")
t_sp<-apc(mod_time_pois, "sp", time_data, type="categorical", nested = c("loc", "sfamily"))
t_pop <-apc(mod_time_pois, "loc", time_data, type = "categorical", nested = c("sp", "sfamily"))
t_t1<-apc(mod_time_pois,"temp1", time_data, type = "binary")
t_t2<-apc(mod_time_pois,"temp2", time_data, type="binary")
t_t3<-apc(mod_time_pois,"temp3", time_data, type="binary")

t_plalan<-apc_withinsp(mod_time_pois, "loc", time_data,type = "categorical", nested="sfamily", sp=4)


#################################
### Making a table #############
################################
apctable<-rbind(cbind(ge_o,ge_s,ge_t1, ge_t2, ge_t3,ge_sp, ge_pop, ge_plalan),
                cbind(t_o,t_s,t_t1,t_t2,t_t3,t_sp, t_pop, t_plalan),
                cbind(gr_o, gr_s, gr_t1,gr_t2,gr_t3,gr_sp, gr_pop, gr_plalan))
apctable<-round(apctable,3)
table1<-as.data.frame(rbind(
  cbind(paste(apctable[1,1],"$\\pm$",apctable[1,2]),
        paste(apctable[1,3], "$\\pm$",apctable[1,4]),
        paste(apctable[1,5], "$\\pm$",apctable[1,6]),
        paste(apctable[1,7], "$\\pm$",apctable[1,8]),                                            
        paste(apctable[1,9], "$\\pm$",apctable[1,10]),
        paste(apctable[1,11], "$\\pm$",apctable[1,12]),
        paste(apctable[1,13], "$\\pm$", apctable[1,14]),
        paste(apctable[1,15], "$\\pm$",apctable[1,16])),
  cbind(paste(apctable[2,1],"$\\pm$",apctable[2,2]),
        paste(apctable[2,3], "$\\pm$",apctable[2,4]),
        paste(apctable[2,5], "$\\pm$",apctable[2,6]),
        paste(apctable[2,7], "$\\pm$",apctable[2,8]),                                            
        paste(apctable[2,9], "$\\pm$",apctable[2,10]),
        paste(apctable[2,11], "$\\pm$",apctable[2,12]), 
        paste(apctable[2,13], "$\\pm$", apctable[2,14]),
        paste(apctable[2,15], "$\\pm$",apctable[2,16])),
  cbind(paste(apctable[3,1],"$\\pm$",apctable[3,2]),
        paste(apctable[3,3], "$\\pm$",apctable[3,4]),
        paste(apctable[3,5], "$\\pm$",apctable[3,6]),
        paste(apctable[3,7], "$\\pm$",apctable[3,8]),                                            
        paste(apctable[3,9], "$\\pm$",apctable[3,10]),
        paste(apctable[3,11], "$\\pm$",apctable[3,12]),
        paste(apctable[3,13], "$\\pm$", apctable[3,14]),
        paste(apctable[3,15], "$\\pm$",apctable[3,16]))
  ))

names(table1)<-c("origin", "stratification", "temperature 1", "temperature 2","temperature 3",
                 "species", "population", "PLALAN population")
row.names(table1)<-c("germination rate (fraction)", "germination date (days)", "growth rate (cm/day)")
table_t<-t(table1)
print(xtable(table_t),tabular.environment = "longtable",include.rownames = TRUE, floating = FALSE, 
      sanitize.text.function = identity)
setwd("C:/Users/Owner/Documents/github/germination_stan")


#############################
######### SAVING ###########
############################

list<-ls()
list<-list[-which(list %in% c("mod_gr","mod_rate","mod_time_pois"))] # list all variables except the one you want to save
save(list = list, file = "apc_germ_balanced.RData")

## ARCHIVE 

# u1.500<-(which(E_u1>250,arr.ind = TRUE))
# E.70<-which(E>80, arr.ind = TRUE)
# u0.500<-(which(E_u0>250, arr.ind = TRUE))
# u1nr<-posterior_predict(mod,newdf1, re.form=NA)
# ppc_dens_overlay(y=time_data$y, yrep = E_u1) +xlim(0,100)
# pdf(file="time_mod_ppchecks.pdf")
# ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois,time_data))+xlim(0,80) + ggtitle("rstanarm")
# ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois_brm_nt,time_data))+xlim(0,80) + ggtitle("brms")
# ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois_brm,time_data))+xlim(0,80) + 
#   ggtitle("brms with truncation")
# dev.off()


# ppc_dens_overlay(y = datax$y,
#                  yrep = E_u1) + xlim(0,100)
# ppc_dens_overlay(y = datax$y,
#                  yrep = E_u0) + xlim(0,100)
# 
# load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois_brm_nt.Rdata")
# load("C:/Users/owner/documents/thesis/stan/mod_time_pois.Rdata")
# pp<-posterior_predict(mod_time_pois_brm)
# pp_s<-as.data.frame(pp[10:20,])
# ppg<-posterior_predict(mod_gr,draws=4000)
# cpp<-colMeans(pp)
# cpg<-colMeans(ppg)
# ppg<-posterior_predict(mod_gr,draws=4000)
# long = melt(pp_s, id.vars= 'var')
# long = melt(as.data.frame(E_u0), id.vars= row.names(E_u0))
# long = melt(as.data.frame(E_u1), id.vars= row.names(E_u1))
# longg = melt(as.data.frame(ppg), id.vars= row.names(ppg))
# 
# 
# ppcheckall<-ggplot(long,aes(value))+geom_density(aes(color=var)) # +  theme(legend.position = "none") +
#   geom_density(aes(rep(datax$y,11)))
# ppcheckall
# ggplot(long[1:120500,], aes(value))+geom_density(aes(color=variable))+xlim(0,50) #+ geom_density(aes(rep(time_data$y,100)))
# plot<-ggplot(long[1:1205000,], aes(value))+geom_density(aes(color=variable),show.legend = FALSE)+xlim(0,70) + 
#   geom_density(aes(rep(time_data$y,1000)),size=2) +
#   geom_density(aes(rep(cpp,1000)), size=2,color="gray")
# plot2<-ggplot(longg[1:1119000,], aes(value))+geom_density(aes(color=variable),show.legend = FALSE)+xlim(0,.5) + 
#   geom_density(aes(rep(datax$y,1000)),size=2) +
#   geom_density(aes(rep(cpg,1000)), size=2,color="gray")
# 
# ggplot(time_data[time_data$origin==0,])+geom_density(aes(y))
# ggplot(long[1:111900,], aes(value))+geom_density(aes(color=variable))+xlim(0,.5) + geom_density(aes(rep(datax$y,100)))
# posterior_vs_prior(mod_gr)
# 
# 
# shinystan::launch_shinystan(mod_time_pois_brm)

