# Fake data for germination date 
# modified from Dan Flynn's "FakeBudburst_Generate_ind.R" file 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two continents, 7 species, two levels of stratification length, and four levels of strating. 2016-04-01 adding interactions. This ends up generating expected differences, but variation in effect sizes across species is minimal currently.
# 2016-05-16 simplifying a lot, but adding individuals. Removed strating and interactions for now.
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Calculating some data parameters: 
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

###### Set up differences for each level
#locdiff = 0 
tempdiff1 = 0.5 #days earlier from 1 to 2
tempdiff2= -1.0 #days earlier from 2 to 3
tempdiff3= 4.5 #days earlier from 3 to 4
origindiff =5
stratdiff = 0.15 #one day later from 30 to 60 

# Generating fake data without interactions first

######## SD for each treatment
tempdiff.sd = 4
origindiff.sd = 4
stratdiff.sd = 1


### Original with interactions below
# interactions. 9 two-way interactions
temporigin1 = -6
temporigin2=-2
temporigin3=-9
tempstrat1 = -0.2
tempstrat2=-0.1
tempstrat3=-0.25
originstrat = -0.15 # 
origintempstrat1 = 0.2
origintempstrat2=0.05
origintempstrat3=0.5


# interactions. 9 two-way interactions
temporigin.sd = 0.01 #
tempstrat.sd = 0.001 # 
originstrat.sd = 0.001 # 
origintempstrat.sd = 0.001

mm <- model.matrix(~(temp+origin+strat)^3, data.frame(temp, origin, strat))
colnames(mm)

#  with individuals

baseinter = 11 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species. 7 species

fake <- vector()

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
    
    coeff <- c(spint[i], 
               rnorm(1, tempdiff1, tempdiff.sd),
               rnorm(1, tempdiff2, tempdiff.sd),
               rnorm(1, tempdiff3, tempdiff.sd),
               rnorm(1, origindiff, origindiff.sd), 
               rnorm(1, stratdiff, stratdiff.sd),
               rnorm(1, temporigin1, temporigin.sd),
               rnorm(1, temporigin2, temporigin.sd),
               rnorm(1, temporigin3, temporigin.sd),
               rnorm(1, tempstrat1, tempstrat.sd),
               rnorm(1, tempstrat2, tempstrat.sd),
               rnorm(1, tempstrat3, tempstrat.sd),
               rnorm(1, originstrat, originstrat.sd),
               rnorm(1, origintempstrat1, origintempstrat.sd),
               rnorm(1, origintempstrat2, origintempstrat.sd),
               rnorm(1, origintempstrat3, origintempstrat.sd)
          
  )
    
    bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
    
    fakex <- data.frame(y=bb, sp = i,
                        temp, origin, strat)
    
    fake <- rbind(fake, fakex)  
  }

summary(lm(y ~ temp*origin*strat, data = fake)) # sanity check 

# save(list=c("fake"), file = "Fake_germdate.RData")
save(fake, file="Fake_germdata.RData")
