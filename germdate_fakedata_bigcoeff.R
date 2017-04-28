## Now with larger coefficients, to show Stan errors 

# Fake data for germination date
# modified from Dan Flynn's "FakeBudburst_Generate_ind.R" file 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two continents, 7 species, two levels of stratification length,  two levels of origin, 
#   and four levels of temperature. Also has interactions and four interactions. 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/undergrads/harold/analyses/germination_stan") 
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

norigin = 2 # origin ==1 is Europe 
ntemp = 4
nstrat = 2

rep = round((nind/(norigin*norigin*ntemp*nsp)), digits=0) # = 11

(ntot = norigin*ntemp*nstrat*rep) #  176 rows

# Build up the data frame
#loc = gl(nloc, rep, length = ntot) #random effect 

temp = gl(ntemp, rep, length = ntot)
origin = as.numeric(as.character(gl(norigin, rep*ntemp, length = ntot, labels=c(0,1))))
strat = as.numeric(as.character(gl(nstrat, rep*ntemp*norigin, length = ntot, c(0,1))))

temp1<-ifelse(temp == 2, 1, 0)
temp2<-ifelse(temp == 3, 1, 0)
temp3<-ifelse(temp == 4, 1, 0)

treatcombo = paste(temp1, temp2, temp3, origin, strat, sep = "_")

(d <- data_frame(temp1, temp2, temp3, origin, strat, treatcombo))

###### Set up differences for each level, for logged response 
#locdiff = 0 
tempdiff1 = 1.5 #logged days earlier from 1 to 2
tempdiff2= -1 #logged days earlier from 2 to 3
tempdiff3= 2 #logged days earlier from 3 to 4
origindiff = 1.9
stratdiff = 1

# Generating fake data without interactions first

######## SD for each treatment
tempdiff.sd = 0.1
origindiff.sd = 0.1
stratdiff.sd = 0.1


### Original with interactions below
# interactions. 10 two-way interactions
temporigin1 = -1
temporigin2 = -2
temporigin3 = -3
tempstrat1 = -2
tempstrat2 = -1
tempstrat3 = -2
originstrat = -3 # 
origintempstrat1 = 3
origintempstrat2 = 3
origintempstrat3 = 3


# interaction sd
temporigin.sd = 0.2 #
tempstrat.sd = 0.1 # 
originstrat.sd = 0.24 # 
origintempstrat.sd = 0.1

mm <- model.matrix(~(temp1+temp2+temp3+origin+strat)^3, data.frame(temp, origin, strat))
mm<-mm[,-grep("temp1:temp2",colnames(mm))]
mm<-mm[,-grep("temp1:temp3", colnames(mm))]
mm<-mm[,-grep("temp2:temp3",colnames(mm))]
colnames(mm)

#  with individuals

baseinter = 2.48 # baseline intercept across all species 
spint <- baseinter + (c(1:nsp)-mean(1:nsp))/2 # different intercepts by species. 7 species

fake <- vector()
cc<- data_frame()

set.seed(100) #so that the created random coefficients can be compared to modeled coefficients 

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
    coeff <- c(spint[i], 
               rnorm(1, tempdiff1, tempdiff.sd),
               rnorm(1, tempdiff2, tempdiff.sd),
               rnorm(1, tempdiff3, tempdiff.sd),
               rnorm(1, origindiff, origindiff.sd), 
               rnorm(1, stratdiff, stratdiff.sd),
               rnorm(1, temporigin1, temporigin.sd),
               rnorm(1, tempstrat1, tempstrat.sd),
               rnorm(1, temporigin2, temporigin.sd),
               rnorm(1, tempstrat2, tempstrat.sd),
               rnorm(1, temporigin3, temporigin.sd),
               rnorm(1, tempstrat3, tempstrat.sd),
               rnorm(1, originstrat, originstrat.sd),
               rnorm(1, origintempstrat1, origintempstrat.sd),
               rnorm(1, origintempstrat2, origintempstrat.sd),
               rnorm(1, origintempstrat3, origintempstrat.sd)
          
  )
    
    bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
    
    fakex <- data.frame(log_y=bb, sp = i,
                        temp1, temp2, temp3, origin, strat)
    
    fake <- rbind(fake, fakex) 
    cc<-rbind(cc, coeff) #to track randomly generated coefficients 
  } 

summary(lm(log_y ~ 
    origin + strat + temp1 + temp2 + temp3 +
    origin*strat + origin*temp1 + origin*temp2 + origin*temp3 +
    strat*temp1 + strat*temp2 + strat*temp3 +
    origin*strat*temp1 +  origin*strat*temp2 + origin*strat*temp3 +
    (1|sp), data = fake)) 


save(fake, file="Fake_germdata_bigcoeff.RData")
