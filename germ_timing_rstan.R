## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())
runstan=FALSE # set to true if running the stan model 

## some good references:
# https://www.r-bloggers.com/bayesian-regression-with-stan-part-1-normal-regression/
# http://m-clark.github.io/docs/IntroBayes.html 

##libraries
library(rstan)
library(shinystan)
library(ggplot2)

##loading files 
setwd("C:/Users/Owner/Documents/GitHub/germination_stan")
load("germs.Rdata") #cleaned and processed data 

#---------------RStan Germination timing ---------------------------------------------

#Setting up the data for the Stan model
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

##add'l random effects -- not in the model yet
# nloc<-length(unique(data$location)) #16 unique locations 
# loc.a<-data$location
# locd<-data.frame(sort(unique(data$location)))
# locd$n<-seq(1,16,1)
# locn<-data.frame()
# for( i in 1:N) {
#   for (j in 1:nloc) {
#     locn[i,j]<-ifelse(loc.a[i]==locd[j,1], j, 0)}
# }
# loc<-data$location
# locn$v17<-rowSums(locn[-1], na.rm=TRUE)
# loc<-locn$v17
# ifelse()
# dummy_variables_fam<-model.matrix(~uniqind, data=data)
# family<-dummy_variables_fam
# nfamily<-length(unique(data$uniqind)) #80 unique seed families 

#putting all the data together: 
datax<-list(N=N, y=y, temp=temp, origin=origin, strat=strat,  nsp=nsp, sp=sp)
            #,nloc=nloc, nfamily=nfamily, loc=loc, family=family)

## fitting the stan model -------------------------------------------------
if (runstan==TRUE) {
#fit <- stan(file = "germdate.stan", data=datax, chains=10, iter=1000) #  divergent transitions above diag 
#fit1 <- stan(file = "germdate.stan", data=datax, chains=10, iter=1000, control = list(adapt_delta = 0.99))
#fit2 <- stan(file = "germdate.stan", data=datax, chains=4, iter=5000, control = list(adapt_delta = 0.99)) #high Rhat, low mixing 
fit3 <- stan(file = "germdate.stan", data=datax, chains=4, iter=20000, warmup=12000, thin=2, 
             control = list(adapt_delta = 0.99)) #This model yields 915 divergent transitions -- all below the diag in the paris plot 
save(fit3, file="germdate_sp_random3.Rdata")
#save(fit2, file = "germdate_sp-random.Rdata")
#save(fit, file="germdate_nore.Rdata")
}

#plotting some posteriors
pdf("fit3_plots.pdf", width=12, height=9)
plot(fit3, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso","mu_b_temp", "mu_b_strat", "mu_b_origin")) 

##trace plots show poor mixing/convergence 
stan_trace(fit3, pars=c("mu_b_inter_to"))
stan_trace(fit3, pars=c("mu_b_inter_ts"))
stan_trace(fit3, pars=c("mu_b_inter_so"))
stan_trace(fit3, pars=c("mu_b_inter_tso"))
stan_trace(fit3, pars=c("mu_b_temp"))
stan_trace(fit3, pars=c("mu_b_strat"))
stan_trace(fit3,pars=c("mu_b_origin")) 

##pair plots show divergent transitions below the diag 

pairs(fit3, pars=c("mu_b_temp", "mu_b_strat", "mu_b_origin"))
pairs(fit3, pars=c("mu_b_inter_to","mu_b_inter_ts", "mu_b_inter_so", "mu_b_inter_tso"))
dev.off()
print(fit3, digits=3)

#----launching shiny stan---------
my_sso <- launch_shinystan(fit3, rstudio = getOption("shinystan.rstudio"))

