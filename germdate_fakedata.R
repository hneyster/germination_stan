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

#rep = 11 # only 1205/(2*2*4*7)=~11 seeds within each combination of treatments. 
rep = round((nind/(norigin*norigin*ntemp)), digits=0) # = 75

(ntot = norigin*ntemp*nstrat*rep) #  176 rows

# Build up the data frame
#loc = gl(nloc, rep, length = ntot) #random effect 

temp = gl(ntemp, rep, length = ntot)
origin = gl(norigin, rep*ntemp, length = ntot)
strat = gl(nstrat, rep*ntemp*norigin, length = ntot)

temp1<-ifelse(temp == 2, 1, 0)
temp2<-ifelse(temp == 3, 1, 0)
temp3<-ifelse(temp == 4, 1, 0)

treatcombo = paste(temp1, temp2, temp3, origin, strat, sep = "_")

(d <- data_frame(temp1, temp2, temp3, origin, strat, treatcombo))

###### Set up differences for each level, for logged response 
#locdiff = 0 
tempdiff1 = 0.01 #logged days earlier from 1 to 2
tempdiff2= -0.01 #logged days earlier from 2 to 3
tempdiff3= 0.02 #logged days earlier from 3 to 4
origindiff = 0.9
stratdiff = 0.01 

# Generating fake data without interactions first

######## SD for each treatment
tempdiff.sd = 0.014
origindiff.sd = 0.4
stratdiff.sd = 0.005


### Original with interactions below
# interactions. 9 two-way interactions
temporigin1 = -0.07
temporigin2 = -0.02
temporigin3 = -0.09
tempstrat1 = -0.002
tempstrat2 = -0.001
tempstrat3 = -0.0025
originstrat = -0.022 # 
origintempstrat1 = 0.005
origintempstrat2 = 0.001
origintempstrat3 = 0.005


# interactions. 9 two-way interactions
temporigin.sd = 0.021 #
tempstrat.sd = 0.0003 # 
originstrat.sd = 0.0004 # 
origintempstrat.sd = 0.001

mm <- model.matrix(~(temp1+temp2+temp3+origin+strat)^3, data.frame(temp, origin, strat))
mm<-mm[,-grep("temp1:temp2",colnames(mm))]
mm<-mm[,-grep("temp1:temp3", colnames(mm))]
mm<-mm[,-grep("temp2:temp3",colnames(mm))]
colnames(mm)

#  with individuals

#setting up random interecpt:
baseinter = 2.48 # baseline intercept across all species 
spint <- baseinter + (c(1:nsp)-mean(1:nsp))/2 # different intercepts by species. 7 species

#setting up random slopes -- not sure I"m doing this right 
sptempdiff1 = seq(-0.005, 0.005, length.out=nsp)
sptempdiff2 = seq(-0.005, 0.005, length.out=nsp)
sptempdiff3 = seq(-0.005, 0.005, length.out=nsp)
sporigindiff =seq(-.6, .6, length.out=nsp) 
spstratdiff = seq(-.01, 0.01, length.out=nsp)
sptemporigin1 = seq(-.005, .005, length.out=nsp)
sptemporigin2 = seq(-.005, .005, length.out=nsp)
sptemporigin3 = seq(-.005, .005, length.out=nsp)
sptempstrat1 = seq(-.00003, .00003, length.out = nsp)
sptempstrat2 = seq(-.00003, .00003, length.out = nsp)
sptempstrat3 = seq(-.00003, .00003, length.out = nsp)
sporiginstrat = seq(-.005, .005, length.out = nsp)
sporigintempstrat1 = seq(-.0001, .0001, length.out = nsp)
sporigintempstrat2 =  seq(-.0001, .0001, length.out = nsp)
sporigintempstrat3 =  seq(-.0001, .0001, length.out = nsp)


fake <- vector()

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
    
    coeff <- c(spint[i], 
               rnorm(1, tempdiff1, tempdiff.sd) +sptempdiff1[i],
               rnorm(1, tempdiff2, tempdiff.sd)+ sptempdiff2[i],
               rnorm(1, tempdiff3, tempdiff.sd) + sptempdiff3[i],
               rnorm(1, origindiff, origindiff.sd) + sporigindiff[i],
               rnorm(1, stratdiff, stratdiff.sd) + spstratdiff[i],
               rnorm(1, temporigin1, temporigin.sd) + sptemporigin1[i],
               rnorm(1, tempstrat1, tempstrat.sd) + sptempstrat1[i], 
               rnorm(1, temporigin2, temporigin.sd) + sptemporigin2[i],
               rnorm(1, tempstrat2, tempstrat.sd) + sptempstrat2[i],
               rnorm(1, temporigin3, temporigin.sd) + sptemporigin3[i],
               rnorm(1, tempstrat3, tempstrat.sd) + sptempstrat3[i],
               rnorm(1, originstrat, originstrat.sd) + sporiginstrat[i],
               rnorm(1, origintempstrat1, origintempstrat.sd) + sporigintempstrat1[i],
               rnorm(1, origintempstrat2, origintempstrat.sd) + sporigintempstrat2[i],
               rnorm(1, origintempstrat3, origintempstrat.sd) + sporigintempstrat3[i]
          
  )
    
    bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
    
    fakex <- data.frame(log_y=bb, sp = i,
                        temp1, temp2, temp3, origin, strat)
    
    fake <- rbind(fake, fakex)  
  }

summary(lm(log_y ~ temp1*temp2*temp3*origin*strat, data = fake)) # sanity check 

# save(list=c("fake"), file = "Fake_germdate.RData")
save(fake, file="Fake_germdata.RData")
