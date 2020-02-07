#Code to plot population-variation for germination rate, time, and growth rate for a single species, PLALAN
#Written by Harold Eyster 1/19/20

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)

setwd("C:/Users/Owner/Documents/GitHub/germination_stan")


##libraries
library(rstan)
library(lme4)
library(ggplot2)
library(rstanarm)
library(tidyr)
library(RCurl)
library(forcats)
library(RColorBrewer)
source("http://peterhaschke.com/Code/multiplot.R") #so that the multiplot function works 
source("https://raw.githubusercontent.com/jaredlander/coefplot/master/R/position.r") # for vertical dodging 

print<-FALSE


####### Germination rate:
#######
#######

#model:
load("C:/Users/Owner/Documents/Thesis/Stan/mod_rate.Rdata") #this is my rstanarm  stanfit object.
m<-mod_rate

# model data: 
load("C:/Users/Owner/Documents/github/germination_stan/rate_data.rdata") 
data<-rate_data
## PLALAN is species 4,
# locations for PLALAN: 
locs<-data[data$sp==4 ,"loc"] %>% unique() %>% sort()
# rate_data %>% group_by(origin) %>% summarise(uniqueloc = n_distinct(sfamily))

# location index: 
load("C:/Users/Owner/Documents/GitHub/germination_stan/germs.Rdata") #cleaned and processed real data
germs.y<-(subset(germs, germinated==1 & sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics
loc<-as.numeric(as.factor(germs.y$loc))
loc_index <- cbind(germs.y$loc, loc) %>% as.data.frame() %>% unique() 


pars_append<-list()
for (i in locs){
  pars <- c(paste("b[origin loc:sp:",i,":4]", sep=""),
            paste("b[strat loc:sp:",i,":4]", sep=""),
            paste("b[temp1 loc:sp:",i,":4]", sep=""),
            paste("b[temp2 loc:sp:",i,":4]", sep=""),
            paste("b[temp3 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp3 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp3 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp3 loc:sp:",i,":4]", sep="")
  )
  pars_append<-(c(pars_append,pars))
}


# This is an Rstanarm object. The coeffs have to be extracted.  There's probably a better way to do this, 
#but I already had this code written up for something else, so used it again.

sum.m <-
  summary(
    m,
    pars = c(
      "b[origin sp:4]",
      "b[strat sp:4]",
      "b[temp1 sp:4]",
      "b[temp2 sp:4]",
      "b[temp3 sp:4]",
      "b[origin:strat sp:4]",
      "b[origin:temp1 sp:4]",
      "b[origin:temp2 sp:4]",
      "b[origin:temp3 sp:4]",
      "b[strat:temp1 sp:4]",
      "b[strat:temp2 sp:4]",
      "b[strat:temp3 sp:4]",
      "b[origin:strat:temp1 sp:4]",
      "b[origin:strat:temp2 sp:4]",
      "b[origin:strat:temp3 sp:4]"
      ,
      pars_append
      
    )
  )

cri.f<-sum.m[,c(1,4,8)] #just selecting the mean and 95% CI.
cri.f_reorder<-cri.f[c(196:210,1:195),] # reording so global plalan is first 
fdf<-data.frame(cri.f_reorder)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rownames(fdf)[1:15], each=length(locs)))), #stdarzing the parameter  names 
    as.numeric(as.character(fdf$mean)),  # the estimate 
    as.numeric(as.character(fdf$X2.5.)), #lower bound, 95% CI
    as.numeric(as.character(fdf$X97.5.)),  #upper bound, 95% CI
    as.numeric(c(rep(1, 15), rep(2, (nrow(fdf)-15)))),  # A variable to signify if the corresponding row is global PLALAN or random loc. 1=global PLALAN, 2=rndm loc
    as.numeric( c( rep(0,15), c(rep(c(locs),15 )))))) #loc variable. Zero when a global 

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "loc") #renaming. 

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Global plalan estimates:
plalan<-c(rep(0, 15), rep(as.numeric(fdf2[c(1:15),"Estimate"]), each=length(locs)))
dff<-fdf2

#adding the fixef estiamtes to the random effect values:
dff$Estimate<-fdf2$Estimate+plalan
dff$`2.5%`<-fdf2$`2.5%`+plalan
dff$`97.5%`<-fdf2$`97.5%`+plalan

dff$var <- fct_inorder(dff$var) #so that the categorical variables plot in the right order 
dff$var2 <- rev(levels(dff$var))
## plotting

pd <- position_dodgev(height = -0.5)
getPalette = colorRampPalette(brewer.pal(13, "Set1"))

locs_names<-loc_index[loc_index$loc %in% locs,]
locs_names<-locs_names[order(as.numeric(as.character(locs_names$loc))),]

fig1<-ggplot(dff, aes(x=Estimate, y=var, color=factor(as.numeric(loc)), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels =c("PLALAN_global", locs_names$V1),
                      values=c("red", getPalette(13)))+ 
  scale_alpha_manual(values=c(1, 0.5))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Germination rate")+
  scale_y_discrete(limits = (unique(sort(dff$var))))+
  theme(legend.position = "none")

#pdf(file.path(figpath, "Fig2.pdf"), width = 7, height = 8)
#fig1
#dev.off()


## Germiantion date 
#model:
load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois.Rdata") #this is my rstanarm  stanfit object.
# model data: 
load("C:/Users/Owner/Documents/github/germination_stan/time_data.rdata") 
## PLALAN is species 4,
# locations for PLALAN: 
locs<-time_data[time_data$sp==4 ,"loc"] %>% unique() %>% sort()

# location index: 
load("C:/Users/Owner/Documents/GitHub/germination_stan/germs.Rdata") #cleaned and processed real data
germs.y<-(subset(germs, germinated==1 & sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics
loc<-as.numeric(as.factor(germs.y$loc))
loc_index <- cbind(data$loc, loc) %>% as.data.frame() %>% unique() 


pars_append<-list()
for (i in locs){
  pars <- c(paste("b[origin loc:sp:",i,":4]", sep=""),
            paste("b[strat loc:sp:",i,":4]", sep=""),
            paste("b[temp1 loc:sp:",i,":4]", sep=""),
            paste("b[temp2 loc:sp:",i,":4]", sep=""),
            paste("b[temp3 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[origin:temp3 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[strat:temp3 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp1 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp2 loc:sp:",i,":4]", sep=""),
            paste("b[origin:strat:temp3 loc:sp:",i,":4]", sep="")
  )
  pars_append<-(c(pars_append,pars))
  print(i)
}


m<-mod_time_pois
# This is an Rstanarm object. The coeffs have to be extracted.  There's probably a better way to do this, 
#but I already had this code written up for something else, so used it again.

sum.m <-
  summary(
    m,
    pars = c(
      "b[origin sp:4]",
      "b[strat sp:4]",
      "b[temp1 sp:4]",
      "b[temp2 sp:4]",
      "b[temp3 sp:4]",
      "b[origin:strat sp:4]",
      "b[origin:temp1 sp:4]",
      "b[origin:temp2 sp:4]",
      "b[origin:temp3 sp:4]",
      "b[strat:temp1 sp:4]",
      "b[strat:temp2 sp:4]",
      "b[strat:temp3 sp:4]",
      "b[origin:strat:temp1 sp:4]",
      "b[origin:strat:temp2 sp:4]",
      "b[origin:strat:temp3 sp:4]"
      ,
      pars_append
      
    )
  )

cri.f<-sum.m[,c(1,4,8)] #just selecting the mean and 95% CI.
cri.f_reorder<-cri.f[c(196:210,1:195),] # reording so global plalan is first 
fdf<-data.frame(cri.f_reorder)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rownames(fdf)[1:15], each=length(locs)))), #stdarzing the parameter  names 
    as.numeric(as.character(fdf$mean)),  # the estimate 
    as.numeric(as.character(fdf$X2.5.)), #lower bound, 95% CI
    as.numeric(as.character(fdf$X97.5.)),  #upper bound, 95% CI
    as.numeric(c(rep(1, 15), rep(2, (nrow(fdf)-15)))),  # A variable to signify if the corresponding row is global PLALAN or random loc. 1=global PLALAN, 2=rndm loc
    as.numeric( c( rep(0,15), c(rep(c(locs),15 )))))) #loc variable. Zero when a global 

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "loc") #renaming. 

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Global plalan  estimates:
plalan<-c(rep(0, 15), rep(as.numeric(fdf2[c(1:15),"Estimate"]), each=length(locs)))
dff<-fdf2

#adding the fixef estiamtes to the random effect values:
dff$Estimate<-fdf2$Estimate+plalan
dff$`2.5%`<-fdf2$`2.5%`+plalan
dff$`97.5%`<-fdf2$`97.5%`+plalan

dff$var <- fct_inorder(dff$var) #so that the categorical variables plot in the right order 
dff$var2 <- rev(levels(dff$var))
## plotting

pd <- position_dodgev(height = -0.5)
getPalette = colorRampPalette(brewer.pal(13, "Set1"))

locs_names<-loc_index[loc_index$loc %in% locs,]
locs_names<-locs_names[order(as.numeric(as.character(locs_names$loc))),]

fig2<-ggplot(dff, aes(x=Estimate, y=var, color=factor(as.numeric(loc)), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels = c("PLALAN_global", locs_names$V1),
                      values=c("red", getPalette(13)))+  scale_alpha_manual(values=c(1, 0.5))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Germination day")+
  scale_y_discrete(limits = (unique(sort(dff$var))))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())




####################
## GROWTH RATE #####
####################
## PLALAN is species 4,
sp_num<-4
#model:
load("C:/Users/Owner/Documents/Thesis/Stan/mod_gr.Rdata") #this is my rstanarm  stanfit object.
m<-mod_gr

# model data: 
load("C:/Users/Owner/Documents/github/germination_stan/datax.Rdata") 
data<-datax
# locations for PLALAN: 
locs<-data[data$sp==sp_num ,"loc"] %>% unique() %>% sort()

# location index: 
load("C:/Users/Owner/Documents/GitHub/germination_stan/germs.Rdata") #cleaned and processed real data
germs.y<-(subset(germs, germinated==1 & sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics
loc<-as.numeric(as.factor(germs.y$loc))
loc_index <- cbind(germs.y$loc, loc) %>% as.data.frame() %>% unique() 


pars_append<-list()
for (i in locs){
  pars <- c(paste("b[origin loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[strat loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[temp1 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[temp2 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[temp3 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:strat loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:temp1 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:temp2 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:temp3 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[strat:temp1 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[strat:temp2 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[strat:temp3 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:strat:temp1 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:strat:temp2 loc:sp:",i,":",sp_num,"]", sep=""),
            paste("b[origin:strat:temp3 loc:sp:",i,":",sp_num,"]", sep="")
  )
  pars_append<-(c(pars_append,pars))
}


# This is an Rstanarm object. The coeffs have to be extracted.  There's probably a better way to do this, 
#but I already had this code written up for something else, so used it again.

sum.m <-
  summary(
    m,
    pars = c(
      paste("b[origin sp:",sp_num,"]",sep=""),
      paste("b[strat sp:",sp_num,"]",sep=""),
      paste("b[temp1 sp:",sp_num,"]",sep=""),
      paste("b[temp2 sp:",sp_num,"]",sep=""),
      paste("b[temp3 sp:",sp_num,"]",sep=""),
      
      paste("b[origin:strat sp:",sp_num,"]",sep=""),
      paste("b[origin:temp1 sp:",sp_num,"]",sep=""),
      paste("b[origin:temp2 sp:",sp_num,"]",sep=""),
      paste("b[origin:temp3 sp:",sp_num,"]",sep=""),
      paste("b[strat:temp1 sp:",sp_num,"]",sep=""),
      paste("b[strat:temp2 sp:",sp_num,"]",sep=""),
      paste("b[strat:temp3 sp:",sp_num,"]",sep=""),
      
      paste("b[origin:strat:temp1 sp:",sp_num,"]",sep=""),
      paste("b[origin:strat:temp2 sp:",sp_num,"]",sep=""),
      paste("b[origin:strat:temp3 sp:",sp_num,"]",sep="")
      ,
      pars_append
      
    )
  )

cri.f<-sum.m[,c(1,4,8)] #just selecting the mean and 95% CI.
cri.f_reorder<-cri.f[c(nrow(cri.f)-14:nrow(cri.f),1:(nrow(cri.f)-15)),] # reording so global plalan is first 
fdf<-data.frame(cri.f_reorder)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rownames(fdf)[1:15], each=length(locs)))), #stdarzing the parameter  names 
    as.numeric(as.character(fdf$mean)),  # the estimate 
    as.numeric(as.character(fdf$X2.5.)), #lower bound, 95% CI
    as.numeric(as.character(fdf$X97.5.)),  #upper bound, 95% CI
    as.numeric(c(rep(1, 15), rep(2, (nrow(fdf)-15)))),  # A variable to signify if the corresponding row is global PLALAN or random loc. 1=global PLALAN, 2=rndm loc
    as.numeric( c( rep(0,15), c(rep(c(locs),15 )))))) #loc variable. Zero when a global 

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "loc") #renaming. 

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Global plalan effect estimates:
plalan<-c(rep(0, 15), rep(as.numeric(fdf2[c(1:15),"Estimate"]), each=length(locs)))
dff<-fdf2

#adding the fixef estiamtes to the random effect values:
dff$Estimate<-fdf2$Estimate+plalan
dff$`2.5%`<-fdf2$`2.5%`+plalan
dff$`97.5%`<-fdf2$`97.5%`+plalan

dff$var <- fct_inorder(dff$var) #so that the categorical variables plot in the right order 

## plotting

pd <- position_dodgev(height = -0.5)
getPalette = colorRampPalette(brewer.pal(13, "Set1"))

locs_names<-loc_index[loc_index$loc %in% locs,]
locs_names<-locs_names[order(as.numeric(as.character(locs_names$loc))),]
locs_names$short <-c("N. France", "Denmark", "S. Slovenia", "S. France", "Austria (high altitude)", "Switzerland",
                     "Germany", "N. Slovenia", "Liechtenstein", "Netherlands", "Austria (low altitude)",
                     "USA--Boston", "USA--Concord")

fig3<-ggplot(dff, aes(x=Estimate, y=var, color=factor(as.numeric(loc)), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual( labels = c("PLALAN_global", locs_names$short),
                      values=c("red", getPalette(13)))+  scale_alpha_manual(values=c(1, 0.5))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Growth rate")+
  scale_y_discrete(limits = (unique(sort(dff$var))))+
  theme(legend.title = element_blank())+
  theme(legend.position = "right", legend.text = element_text(size = 8))+
  theme(axis.text.y = element_blank())

if(print==TRUE){
  png(filename = "PLALAN_pops_plot.png", width = 14, height= 10, units = "in", res = 600)
  multiplot(fig1,fig2,fig3, cols=3)
  dev.off()
}


### Climate data ### 
clim_raw<-read.csv("avgtemps.csv")
clim<-merge(clim_raw, locs_names, by.x="loc_name", by.y="V1", all = FALSE) ## just the locs with PLALAN in them
names(clim)[7] <- "loc_num"
pd2<- -.2
 
fig4<-  ggplot(data=clim) + 
    geom_point(position=position_dodge(-.2), aes(y=clim$Apr, x=1, color=factor(as.numeric(clim$loc_num))), size=8, shape=15)+
    geom_point(position=position_dodge(-.2),aes(x=2, y=clim$Apr, color=factor(as.numeric(clim$loc_num))), size=8, shape=16)+
    geom_point(position=position_dodge(-.3),aes(x=3, y=clim$May, color=factor(as.numeric(clim$loc_num))), size=8, shape=17)+
    scale_color_manual(labels=clim$short, values=getPalette(13))+
    scale_alpha_manual(values=c(.8))+
    scale_x_continuous(breaks=c(1,2,3), labels=c("Mar","Apr","May")) + 
    ylab("Mean temperature (C)")+
    xlab("Month")+
  theme(legend.title = element_blank()) 

## NOw all the locs: 
al<-.8
clim_all<-merge(clim_raw, loc_index, by.x="loc_name", by.y="V1", all=FALSE) 
clim_all$short <-c("Norway", locs_names$short[1:9], "W. France",
                   locs_names$short[10:12], "USA--Harvard Forest", locs_names$short[13])
names(clim_all)[7] <- "loc_num"
clim<-clim_all
fig5<-  ggplot(data=clim) + 
  geom_point(position=position_dodge(-.2), aes(y=clim$Apr, x=1, color=factor(as.numeric(clim$loc_num))), size=8, shape=16, alpha=al)+
  geom_point(position=position_dodge(-.2),aes(x=2, y=clim$Apr, color=factor(as.numeric(clim$loc_num))), size=8, shape=16, alpha=al)+
  geom_point(position=position_dodge(-.3),aes(x=3, y=clim$May, color=factor(as.numeric(clim$loc_num))), size=8, shape=16, alpha=al)+
  scale_color_manual(labels=clim$short, values=getPalette(nrow(clim)))+
  scale_x_continuous(breaks=c(1,2,3), labels=c("Mar","Apr","May")) + 
  ylab(expression(paste("Mean temperature 1970-2000 (",~degree*C,")")))+
  xlab("Month")+
  theme(legend.title = element_blank()) 

al<-.5
fig6<-  ggplot(data=clim) + 
  geom_point(position=position_dodge(-.2), aes(y=clim$Apr, x=1, color=factor(as.numeric(clim$loc_num))), size=8, shape=16,alpha=al)+
  geom_point(position=position_dodge(-.2),aes(x=2, y=clim$Apr, color=factor(as.numeric(clim$loc_num))), size=8, shape=16,alpha=al)+
  geom_point(position=position_dodge(-.3),aes(x=3, y=clim$May, color=factor(as.numeric(clim$loc_num))), size=8, shape=16, alpha=al)+
  geom_point(position=position_dodge(-.2), aes(x=4, y=(clim$May-clim$Apr),
     color=factor(as.numeric(clim$loc_num))), size=8,shape=17, alpha=.3)+
  scale_color_manual(labels=clim$short, values=c(rep("red",13), rep("blue",3)))+
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("Mar","Apr","May", "May-Apr")) + 
  ylab("Mean temperature (C)")+
  xlab("Month")+
  theme(legend.title = element_blank()) 
#   guides(color = guide_legend(override.aes = list(size=5)))

# sum.m <-
#   summary(
#     m,
#     pars = c(
#       "(Intercept)",
#       "origin",
#       "strat",
#       "temp1",
#       "temp2",
#       "temp3",
#       "origin:strat",
#       "origin:temp1",
#       "origin:temp2",
#       "origin:temp3",
#       "strat:temp1",
#       "strat:temp2",
#       "strat:temp3",
#       "origin:strat:temp1",
#       "origin:strat:temp2",
#       "origin:strat:temp3"
#       
#       ,
#       "b[origin sp:4]",
#       "b[strat sp:4]",
#       "b[temp1 sp:4]",
#       "b[temp2 sp:4]",
#       "b[temp3 sp:4]",
#       "b[origin:strat sp:4]",
#       "b[origin:temp1 sp:4]",
#       "b[origin:temp2 sp:4]",
#       "b[origin:temp3 sp:4]",
#       "b[strat:temp1 sp:4]",
#       "b[strat:temp2 sp:4]",
#       "b[strat:temp3 sp:4]",
#       "b[origin:strat:temp1 sp:4]",
#       "b[origin:strat:temp2 sp:4]",
#       "b[origin:strat:temp3 sp:4]"
#       ,
#       pars_append
#       
#     )
#   )
# 
# cri.f<-sum.m[c(2:nrow(sum.m)),c(1,4,8)] #just selecting the mean and 95% CI. Removing the intercept 
# fdf<-data.frame(cri.f)
# #binding 
# fdf2<-as.data.frame(
#   cbind(
#     (c(rownames(fdf)[1:15],rownames(fdf)[1:15], rep(rev(rownames(fdf)[1:15]), each=length(locs)))), #stdarzing the parameter  names 
#     as.numeric(as.character(fdf$mean)),  # the estimate 
#     as.numeric(as.character(fdf$X2.5.)), #lower bound, 95% CI
#     as.numeric(as.character(fdf$X97.5.)),  #upper bound, 95% CI
#     as.numeric(c(rep(1, 15), rep(2, 15), rep(3, (nrow(fdf)-30)))),  # A variable to signify if the corresponding row is a fixed  or random effect. 1=global, 2=rndm sp, 3 = rndm loc
#     as.numeric( c(rep(0,15), c(rep(1,15)), c(rep(seq(2:(length(locs)+1)),15 )))))) #loc variable. Zero when a global, 1 when for sp mean 

