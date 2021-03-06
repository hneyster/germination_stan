## Average Predictive Comparisons for plant germiantion 
# updated 1/20/20 to enable it to compare populations 
options(stringsAsFactors = FALSE)
library(dplyr)
library(arm)
library(rstan)
library(rstanarm)
library(reshape2)

################################
####### The function ###########
################################

apc<-function(mod,u,df,type = c("binary","numerical", "categorical"), nested, freq_in_data=FALSE){ 
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
    #sample_draws<-sample(1:S,1000)
    #E_u1<-E_u1[sample_draws,] 
    #E_u0<-E_u0[sample_draws,]
    E_diff<-((E_u1-E_u0)^2) #the squared difference, as in equation 6. 
    sum_theta<-colSums(E_diff) #summing across model draws 
    sum_nu_theta<-sum(sum_theta) #summing across nu
    num<-sum_nu_theta #this is the numerator in equation 6 
    denom<-S*n #this is the denominator in equation 6. 
    APC<-(num/denom)^(1/2) #this is the average predictive comparison 
    
    #now calculating the stadard error
    sum_nu<-rowSums(E_diff)/n
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
      #categorical rndm effects may be  correlated with nested rndm effects
      #Here, we retain these relations: 
      subdf<-data.frame((newdf[newdf[u]==i,nested])) #identifying which nested effects correlate
      if (freq_in_data==FALSE){
        subdf<- unique(subdf) 
        }# should we create the sample with the same 
      # data frequency as the data (freq_in_data=TRUE)? Or with balanced frequency (freq_in_data=FALSE)? 
      newdf1[((n*(i-1))+1):(n*i),nested] <- subdf[sample(nrow(subdf),n,replace=TRUE), ] 
      #then sampling from these with replacement, and inserting back in the alternative dataframe
    }
    E_u1<-posterior_predict(mod,newdata=newdf1,draws=S,seed=248) #each of these columns
    #represents the expected value for a different value of nu; each set of n columns represents the expected values
    # for each value of upsilon. each row represents the expected value according to a different model draw
    E_u0<-posterior_predict(mod,newdata=newdf0,draws=S,seed=248) # Now the same, but for the unadulterated data 
    sample_draws<-sample(1:S,1000)
    E_u1<-E_u1[sample_draws,] 
    E_u0<-E_u0[sample_draws,]
    E_diff<-((E_u1-E_u0)^2)
    sum_theta<-colSums(E_diff) #summing across model draws 
    sum_nu_u_theta<-sum(sum_theta)
    num<-sum_nu_u_theta
    denom<-S*n*K
    APC<-(num/denom)^(1/2)
    
    #now calculating SE: 
    sum_nu_u<-rowSums(E_diff)/(n*K)
    apc_vec<-rep(APC^2,S) # a vector of apc 
    SE<-(1/(2*APC))*(sqrt(1/(S-1)*(sum  ((sum_nu_u-apc_vec)^2)  ))) #The standard error
    #return(paste('APC =',APC, ', SE = ',SE))
    return(data.frame("APC"=APC, "SE"=SE))
  }
  else {
    return(paste('variable of interest must be binary or categorical'))
  }
}

