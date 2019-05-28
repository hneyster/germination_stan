rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(arm)
library(rstan)
library(rstanarm)
library(reshape2)
library(ggplot2)
library(xtable)


################################
####### The function ###########
################################

apc<-function(mod,u,df,type = c("binary","numerical", "categorical")){ 
  S<- ifelse((attributes(mod)$class[1]=="brmsfit"), 4000, 1000) #brms objects predict from the entire suite of model 
  #draws (4000 in this case), so we'll have to use S = 4000 for the brms objects 
  #mod= model, u= the parameter of interest, df=dataframe, type = type of input variable of interest 
  n<-nrow(df) 
  newdf<-df[,c("temp1","temp2", "temp3","origin","strat","sp","loc","sfamily")] #subsetting just the observation data, nu and upsilon 
  if (type=="binary") {
    newdf1<-newdf
    newdf1[[u]]<-rep(1,n)
    newdf0<-newdf
    newdf0[[u]]<-rep(0,n)
    E_u1<-posterior_predict(mod,newdata=newdf1,draws=S,seed=248) #each of these columns
    #represents the expected value for a different value of nu; each row represents the 
    #expected value according to a different model draw; each of these is for upsilon=1 
    E_u0<-posterior_predict(mod,newdata=newdf0,draws=S,seed=248) # Now the same, but for upsilon=0
    sample_draws<-sample(1:S,1000)
    E_u1<-E_u1[sample_draws,] 
    E_u0<-E_u0[sample_draws,]
    E_u1<-ifelse(E_u1<134, E_u1,NA) 
    E_u0<-ifelse(E_u0<134, E_u0,NA)
    E_diff<-((E_u1-E_u0)^2) #the squared difference, as in equation 6. 
    sum_theta<-colSums(E_diff, na.rm=TRUE) #summing across model draws 
    sum_nu_theta<-sum(sum_theta,na.rm = TRUE) #summing across nu
    num<-sum_nu_theta #this is the numerator in equation 6 
    denom<-S*n #this is the denominator in equation 6. 
    APC<-(num/denom)^(1/2) #this is the average predictive comparison 
    
    #now calculating the stadard error
    sum_nu<-rowSums(E_diff,na.rm = TRUE)/n
    apc_vec<-rep(APC^2,S) # a vector of apc 
    SE<-(1/(2*APC))*(sqrt(1/(S-1)*(sum  ((sum_nu-apc_vec)^2)  ))) #The standard error
    #return(paste('APC =',APC, ', SE = ',SE))
    return(data.frame("APC"=APC, "SE"=SE))
    
  }
  else if (type=="categorical"){
    K=length(unique(df[[u]])) #the number of values that u takes 
    u_k<-rep(seq(1,K),each=n) #This is a vector of values of u, so that each value of observation has every value of u. 
    #has dimensions n*K
    newdf0<-newdf[rep(row.names(newdf), K),1:ncol(newdf)] #this takes the data and replicates it one time for each value of u
    newdf1<-newdf0 #this is just the replicated data 
    newdf1[[u]]<-u_k #inserting the new values of upsilon 
    for (i in 1:K) {
      subdf<-(newdf[newdf$sp==i,7:8])
      newdf1[((n*(i-1))+1):(n*i),7:8] <- subdf[sample(nrow(subdf),n,replace=TRUE), ]
    }
    E_u1<-posterior_predict(mod,newdata=newdf1,draws=S,seed=248) #each of these columns
    #represents the expected value for a different value of nu; each set of n columns represents the expected values
    # for each value of upsilon. each row represents the expected value according to a different model draw
    E_u0<-posterior_predict(mod,newdata=newdf0,draws=S,seed=248) # Now the same, but for the unadulterated data 
    sample_draws<-sample(1:S,1000)
    E_u1<-E_u1[sample_draws,] 
    E_u0<-E_u0[sample_draws,]
    E_u1<-ifelse(E_u1<134, E_u1,NA) 
    E_u0<-ifelse(E_u0<134, E_u0,NA)
    E_diff<-((E_u1-E_u0)^2)
    sum_theta<-colSums(E_diff, na.rm=TRUE) #summing across model draws 
    sum_nu_u_theta<-sum(sum_theta, na.rm=TRUE)
    num<-sum_nu_u_theta
    denom<-S*n*K
    APC<-(num/denom)^(1/2)
    
    #now calculating SE: 
    sum_nu_u<-rowSums(E_diff,na.rm=TRUE)/(n*K)
    apc_vec<-rep(APC^2,S) # a vector of apc 
    SE<-(1/(2*APC))*(sqrt(1/(S-1)*(sum  ((sum_nu_u-apc_vec)^2)  ))) #The standard error
    #return(paste('APC =',APC, ', SE = ',SE))
    return(data.frame("APC"=APC, "SE"=SE))
  }
  else {
    return(paste('variable of interest must be binary or categorical'))
  }
}


##################################
#### APCs for Growth Rate #### 
#################################
load("C:/Users/Owner/Documents/Thesis/Stan/mod_gr.Rdata")
load("C:/Users/Owner/Documents/github/germination_stan/datax.Rdata") 
gr_o<-apc(mod_gr,"origin",datax,type = "binary")
gr_s<-apc(mod_gr,"strat",datax,type="binary")
gr_sp<-apc(mod_gr, "sp",datax,type="categorical")

###################################
### APCs for germ rate model######## 
###################################


load("C:/Users/Owner/Documents/Thesis/Stan/mod_rate.Rdata") 
load("C:/Users/Owner/Documents/github/germination_stan/rate_data.Rdata") 
ge_o<-apc(mod_rate,"origin",rate_data,type="binary")
ge_s<-apc(mod_rate,"strat",rate_data,type="binary")
ge_sp<-apc(mod_rate,"sp",rate_data,type="categorical")

###################################
#### APCs for germ timing model ###
###################################

load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois_brm.Rdata") #this is a brms stanfit object. 
load("C:/Users/Owner/Documents/github/germination_stan/time_data.rdata") 
t_o<-apc(mod_time_pois_brm,"origin",time_data,type = "binary")
t_s<-apc(mod_time_pois,"strat",time_data,type = "binary")
t_sp<-apc(mod_time_pois_brm, "sp", time_data, type="categorical")

#################################
### Making a table #############
################################
apctable<-rbind(cbind(ge_o,ge_s,ge_sp),cbind(t_o,t_s,t_sp),cbind(gr_o, gr_s, gr_sp))
apctable<-round(apctable,3)
table1<-as.data.frame(rbind(
cbind(paste(apctable[1,1],"$\\pm$",apctable[1,2]),
      paste(apctable[1,3], "$\\pm$",apctable[1,4]),
      paste(apctable[1,5], "$\\pm$",apctable[1,6])), 
cbind(paste(apctable[2,1],"$\\pm$",apctable[2,2]),
      paste(apctable[2,3], "$\\pm$",apctable[2,4]),
      paste(apctable[2,5], "$\\pm$",apctable[2,6])),  
cbind(paste(apctable[3,1],"$\\pm$",apctable[3,2]),
      paste(apctable[3,3], "$\\pm$",apctable[3,4]),
      paste(apctable[3,5], "$\\pm$",apctable[3,6]))    
))
names(table1)<-c("origin", "stratification", "species")
row.names(table1)<-c("germination rate", "germination date", "growth rate")
print(xtable(table1),tabular.environment = "longtable",include.rownames = TRUE, floating = FALSE, 
      sanitize.text.function = identity)

## ARCHIVE 

u1.500<-(which(E_u1>250,arr.ind = TRUE))
E.70<-which(E>80, arr.ind = TRUE)
u0.500<-(which(E_u0>250, arr.ind = TRUE))
u1nr<-posterior_predict(mod,newdf1, re.form=NA)
ppc_dens_overlay(y=time_data$y, yrep = E_u1) +xlim(0,100)
pdf(file="time_mod_ppchecks.pdf")
ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois,time_data))+xlim(0,80) + ggtitle("rstanarm")
ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois_brm_nt,time_data))+xlim(0,80) + ggtitle("brms")
ppc_dens_overlay(y=time_data$y, yrep=posterior_predict(mod_time_pois_brm,time_data))+xlim(0,80) + 
  ggtitle("brms with truncation")
dev.off()


# ppc_dens_overlay(y = datax$y,
#                  yrep = E_u1) + xlim(0,100)
# ppc_dens_overlay(y = datax$y,
#                  yrep = E_u0) + xlim(0,100)
# 
 load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois_brm_nt.Rdata")
 load("C:/Users/owner/documents/thesis/stan/mod_time_pois.Rdata")
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

