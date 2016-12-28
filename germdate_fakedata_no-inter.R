# Fake data for germination date 
# modified from Dan Flynn's "FakeBudburst_Generate_ind.R" file 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two continents, 7 species, two levels of stratification length,  two levels of origin, 
#   and four levels of temperature. Also has interactions and four interactions. 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

if(length(grep("danflynn", getwd())>0)) { 
  setwd("~/Documents/git/germination_stan") 
} else 
  setwd("C:/Users/Owner/Documents/GitHub/germination_stan")

# germidn<-subset(germid, sp!="PLAMED" & sp!="PLACOR")
# length(unique(germidn$uniqueid))
# mean(subset(germidn, germinated==1 & origin=="USA")$daysfromstart)-mean(subset(germidn, germinated==1 & origin=="Europe")$daysfromstart)

library("dplyr")

nloc = 16
nsp = 7
nfamily = 80
nind=1205

norigin=2 # origin ==1 is Europe 
ntemp = 4
nstrat = 2

rep = 11 # only 1205/(2*2*4*7)=~11 seeds within each combination of treatments. 

(ntot = norigin*ntemp*nstrat*rep) # 1200 rows

# Build up the data frame
#loc = gl(nloc, rep, length = ntot) #random effect 

temp = gl(ntemp, rep, length = ntot)
origin = gl(norigin, rep*ntemp, length = ntot)
strat = gl(nstrat, rep*ntemp*norigin, length = ntot)

treatcombo = paste(temp, origin, strat, sep = "_")

(d <- data_frame(temp, origin, strat, treatcombo))

###### Set up differences for each level, for logged response 
#locdiff = 0 
tempdiff1 = 0.01 #logged days earlier from 1 to 2
tempdiff2= -0.01 #logged days earlier from 2 to 3
tempdiff3= .02 #logged days earlier from 3 to 4
origindiff =0.9
stratdiff = 0.01 

# Generating fake data without interactions first

######## SD for each treatment
tempdiff.sd = 0.014
origindiff.sd = 0.4
stratdiff.sd = .005


mm <- model.matrix(~(temp+origin+strat), data.frame(temp, origin, strat))
colnames(mm)

#  with individuals

baseinter = 2.48 # baseline intercept across all species 
spint <- baseinter + (c(1:nsp)-mean(1:nsp))/2 # different intercepts by species. 7 species

fake <- vector()

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
    
    coeff <- c(spint[i], 
               rnorm(1, tempdiff1, tempdiff.sd),
               rnorm(1, tempdiff2, tempdiff.sd),
               rnorm(1, tempdiff3, tempdiff.sd),
               rnorm(1, origindiff, origindiff.sd), 
               rnorm(1, stratdiff, stratdiff.sd))
               
          
  
    
    bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
    
    fakex <- data.frame(log_y=bb, sp = i,
                        temp, origin, strat)
    
    fake <- rbind(fake, fakex)  
  }

summary(lm(log_y ~ temp+origin+strat, data = fake)) # sanity check 

# save(list=c("fake"), file = "Fake_germdate.RData")
save(fake, file="Fake_germdata_no-inter.RData")
