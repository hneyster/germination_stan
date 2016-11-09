# Fake data for germination date 
# modified from Dan Flynn's "FakeBudburst_Generate_ind.R" file 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two continents, 7 species, two levels of stratification length, and four levels of strating. 2016-04-01 adding interactions. This ends up generating expected differences, but variation in effect sizes across species is minimal currently.
# 2016-05-16 simplifying a lot, but adding individuals. Removed strating and interactions for now.
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Calculating some data parameters: 
if(length(grep("danflynn", getwd())>0)) { 
  setwd("~/Documents/git/eysterthesis") 
} else 
  setwd("C:/Users/Owner/Documents/GitHub/eysterthesis")

# germidn<-subset(germid, sp!="PLAMED" & sp!="PLACOR")
# length(unique(germidn$uniqueid))
# mean(subset(germidn, germinated==1 & origin=="USA")$daysfromstart)-mean(subset(germidn, germinated==1 & origin=="Europe")$daysfromstart)

# install.packages("dplyr")
library("dplyr")

nsite = 27
nsp = 7
npop= 2
nind = 3401

norigin=2 # origin ==1 is Europe 
ntemp = 4
nstrat = 2

rep = 1 # only 1 individual within each combination of treatments. 

(ntot = nsite*ntemp*nstrat*rep) # 216 rows

# Build up the data frame
site = gl(nsite, rep, length = ntot)

temp = gl(ntemp, rep*nsite, length = ntot)
origin = gl(norigin, rep*nsite*ntemp, length = ntot)
strat = gl(nstrat, rep*nsite*ntemp, length = ntot)

treatcombo = paste(temp, origin, strat, sep = "_")

(d <- data_frame(site, temp, origin, strat, treatcombo))

###### Set up differences for each level
sitediff = 0 
tempdiff = -5 # days earlier from 1 to 2
origindiff = 3.5
stratdiff = -3.5

# Generating fake data without interactions first

######## SD for each treatment
sitediff.sd = 3 
tempdiff.sd = 1 
origindiff.sd = 1
stratdiff.sd = 1.5

mm <- model.matrix(~site+temp+origin+strat, data.frame(site, temp, origin))
colnames(mm) # *** Look at the colnames here. This is the order you will have to provide the 'coeff' for in the loop below. See that there are multiple columns for the treatments with multiple levels

#  with individuals

baseinter = 10.4 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species. 7 species

fake <- vector()

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
  
  # this will have to be spint[i] still as the intercept, but then coeffs for site2, site3... site 27, temp2, temp3, temp4, origindiff, stratdiff. The two-level treatments are ok, just site and temp need to be expanded.
  
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             
             rnorm(1, tempdiff, tempdiff.sd),
             rnorm(1, origindiff, origindiff.sd), 
             rnorm(1, stratdiff, stratdiff.sd)
             ##
  )
  
  bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, ind = paste(i, j, sep="_"),
                      site, temp, origin, strat)
  
  fake <- rbind(fake, fakex)  
}
} #getting an error: "Error in mm %*% coeff : non-conformable arguments" 

summary(lm(bb ~ (site+temp+origin+strat)^2, data = fake)) # sanity check 

save(list=c("fake"), file = "Fake_germdate.RData")



### Original with interactions below
# interactions. 9 two-way interactions
sitetemp = 0
siteorigin = 0
sitestrat = 0
temporigin = 3.5 # positive 3.5. So at the temp level, the effect of origin is muted by 3.5 days.
tempstrat = 11 # both positive ~ 10. 
originstrat = 0.1 # 

######## SD for each treatment
sitediff.sd = 3 
tempdiff.sd = 1 
origindiff.sd = 1
stratdiff.sd = 1.5


# interactions. 9 two-way interactions
sitetemp.sd = 1
siteorigin.sd = 1 
sitestrat.sd = 2
temporigin.sd = 1
tempstrat.sd = 1.5
originstrat.sd = 1

mm <- model.matrix(~(site+temp+origin+strat)^2, data.frame(site, temp, origin))
colnames(mm)

#  with individuals

baseinter = 10.4 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species. 7 species

fake <- vector()

for(i in 1:nsp){ # loop over species. i = 1
  
  # Give species different difference values, drawn from normal.
    
    coeff <- c(spint[i], 
               rnorm(1, sitediff, sitediff.sd),
               rnorm(1, tempdiff, tempdiff.sd),
               rnorm(1, origindiff, origindiff.sd), 
               rnorm(1, stratdiff, stratdiff.sd),
               rnorm(1, sitetemp, sitetemp.sd), 
               rnorm(1, siteorigin, siteorigin.sd),
               rnorm(1, sitestrat, sitestrat.sd),
               rnorm(1, temporigin, temporigin.sd),
               rnorm(1, tempstrat, tempstrat.sd),
               rnorm(1, originstrat, originstrat.sd)
  )
    
    bb <- rnorm(n = length(temp), mean = mm %*% coeff, sd = 0.1)
    
    fakex <- data.frame(bb, sp = i, ind = paste(i, j, sep="_"),
                        site, temp, origin, strat)
    
    fake <- rbind(fake, fakex)  
  }
} #getting an error: "Error in mm %*% coeff : non-conformable arguments" 

summary(lm(bb ~ (site+temp+origin+strat)^2, data = fake)) # sanity check 

save(list=c("fake"), file = "Fake_germdate.RData")
