#Code to plot model coefficients for germination rate, time, and growth rate. 
#Written by Harold Eyster 1/15/18
# revised with flipped axes (variables on y-axis) 7/22/20

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())

source("http://peterhaschke.com/Code/multiplot.R") #so that the multiplot function works 
#source("https://raw.githubusercontent.com/jaredlander/coefplot/master/R/position.r") # for vertical dodging 



##libraries
library(rstan)
library(lme4)
library(ggplot2)
library(rstanarm)
library(brms)
library(tidyr)
library(RCurl)
library(forcats)
library(RColorBrewer)
library(here)
library(cowplot)

# colors: 
getPalette = colorRampPalette(brewer.pal(7, "Set1"))

# plotting function: 
modplot<-function(model,title, units, legend){
  

# This is an Rstanarm object. The coeffs have to be extracted.  There's probably a better way to do this, 
#but I already had this code written up for something else, so used it again.
sum.m <-
  summary(probs=c(0.025,0.975),
    model,
    pars = c(
      "(Intercept)",
      "origin",
      "strat",
      "temp1",
      "temp2",
      "temp3",
      "origin:strat",
      "origin:temp1",
      "origin:temp2",
      "origin:temp3",
      "strat:temp1",
      "strat:temp2",
      "strat:temp3",
      "origin:strat:temp1",
      "origin:strat:temp2",
      "origin:strat:temp3"
      ,
      "b[origin sp:1]",
      "b[strat sp:1]",
      "b[temp1 sp:1]",
      "b[temp2 sp:1]",
      "b[temp3 sp:1]",
      "b[origin:strat sp:1]",
      "b[origin:temp1 sp:1]",
      "b[origin:temp2 sp:1]",
      "b[origin:temp3 sp:1]",
      "b[strat:temp1 sp:1]",
      "b[strat:temp2 sp:1]",
      "b[strat:temp3 sp:1]",
      "b[origin:strat:temp1 sp:1]",
      "b[origin:strat:temp2 sp:1]",
      "b[origin:strat:temp3 sp:1]"
      ,
      "b[origin sp:2]",
      "b[strat sp:2]",
      "b[temp1 sp:2]",
      "b[temp2 sp:2]",
      "b[temp3 sp:2]",
      "b[origin:strat sp:2]",
      "b[origin:temp1 sp:2]",
      "b[origin:temp2 sp:2]",
      "b[origin:temp3 sp:2]",
      "b[strat:temp1 sp:2]",
      "b[strat:temp2 sp:2]",
      "b[strat:temp3 sp:2]",
      "b[origin:strat:temp1 sp:2]",
      "b[origin:strat:temp2 sp:2]",
      "b[origin:strat:temp3 sp:2]"
      ,
      "b[origin sp:3]",
      "b[strat sp:3]",
      "b[temp1 sp:3]",
      "b[temp2 sp:3]",
      "b[temp3 sp:3]",
      "b[origin:strat sp:3]",
      "b[origin:temp1 sp:3]",
      "b[origin:temp2 sp:3]",
      "b[origin:temp3 sp:3]",
      "b[strat:temp1 sp:3]",
      "b[strat:temp2 sp:3]",
      "b[strat:temp3 sp:3]",
      "b[origin:strat:temp1 sp:3]",
      "b[origin:strat:temp2 sp:3]",
      "b[origin:strat:temp3 sp:3]"
      ,
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
      "b[origin sp:5]",
      "b[strat sp:5]",
      "b[temp1 sp:5]",
      "b[temp2 sp:5]",
      "b[temp3 sp:5]",
      "b[origin:strat sp:5]",
      "b[origin:temp1 sp:5]",
      "b[origin:temp2 sp:5]",
      "b[origin:temp3 sp:5]",
      "b[strat:temp1 sp:5]",
      "b[strat:temp2 sp:5]",
      "b[strat:temp3 sp:5]",
      "b[origin:strat:temp1 sp:5]",
      "b[origin:strat:temp2 sp:5]",
      "b[origin:strat:temp3 sp:5]"
      ,
      "b[origin sp:6]",
      "b[strat sp:6]",
      "b[temp1 sp:6]",
      "b[temp2 sp:6]",
      "b[temp3 sp:6]",
      "b[origin:strat sp:6]",
      "b[origin:temp1 sp:6]",
      "b[origin:temp2 sp:6]",
      "b[origin:temp3 sp:6]",
      "b[strat:temp1 sp:6]",
      "b[strat:temp2 sp:6]",
      "b[strat:temp3 sp:6]",
      "b[origin:strat:temp1 sp:6]",
      "b[origin:strat:temp2 sp:6]",
      "b[origin:strat:temp3 sp:6]"
      ,
      "b[origin sp:7]",
      "b[strat sp:7]",
      "b[temp1 sp:7]",
      "b[temp2 sp:7]",
      "b[temp3 sp:7]",
      "b[origin:strat sp:7]",
      "b[origin:temp1 sp:7]",
      "b[origin:temp2 sp:7]",
      "b[origin:temp3 sp:7]",
      "b[strat:temp1 sp:7]",
      "b[strat:temp2 sp:7]",
      "b[strat:temp3 sp:7]",
      "b[origin:strat:temp1 sp:7]",
      "b[origin:strat:temp2 sp:7]",
      "b[origin:strat:temp3 sp:7]"
    )
  )

cri.f<-sum.m[c(2:121),c(1,4,5)] #just selecting the mean and 95% CI. Removing the intercept 
fdf<-data.frame(cri.f)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rev(rownames(fdf)[1:15]), each=7))), #stdarzing the parameter  names 
    as.numeric(as.character(fdf$mean)),  # the estimate 
    as.numeric(as.character(fdf$X2.5.)), #lower bound, 95% CI
    as.numeric(as.character(fdf$X97.5.)),  #upper bound, 95% CI
    as.numeric(c(rep(1, 15), rep(2, 105))),  # A variable to signify if the corresponding row is a fixed  or random effect. 1=global, 2=rndm
    as.numeric( c(rep(0,15), rep(seq(1,7),15 ))))) #sp variable. Zero when a factor 

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "sp") #renaming. 

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Fixed effect estimates:
fixed<-c(rep(0, 15), rep(as.numeric(rev(fdf2[c(1:15),2])), each=7))
dff<-fdf2

#adding the fixef estiamtes to the random effect values:
dff$Estimate<-fdf2$Estimate+fixed
dff$`2.5%`<-fdf2$`2.5%`+fixed
dff$`97.5%`<-fdf2$`97.5%`+fixed

dff$var <- fct_inorder(dff$var)# %>% fct_rev() # so that the categorical variables plot in the right order 
dff$rndm <- fct_rev(dff$rndm) 
dff$sp <- factor(dff$sp, levels=c(0,1:3,4:7))

dff<-dff[c(16:120, 1:15),] # so that the main effects plot on top 
pd <- position_dodge(width =  0.5)

var_names <-c("origin","strat","temp1" ,"temp2","temp3","origin × strat","origin × temp1","origin × temp2","origin × temp3",
              "strat × temp1", "strat × temp2","strat × temp3","origin × strat × temp1", "origin × strat × temp2",
              "origin × strat × temp3")
fig<-
    
    ggplot(dff, aes(y=Estimate, x=var, color=(sp),alpha=rndm, size=factor(rndm))) + # , alpha=factor(r,ndm)))+
    geom_errorbar(aes(ymin=(`2.5%`), ymax=(`97.5%`)), position=pd, size=.8, width =0)+
    geom_point(position=pd)+
    #geom_point(data = subset(dff, rndm == 1), position=pd)+
    geom_hline(yintercept=0)+
    scale_colour_manual(labels = c("Global mean", expression(italic("Capsella bursa-pastoris")), 
                                   expression(italic("Chelidonium majus")), expression(italic("Dactylis glomerata")),
                                   expression(italic("Plantago lanceolata")), expression(italic("Plantago major")),
                                   expression(italic("Rumex crispus")),expression(italic("Taraxacum officinale"))),
                        values=c("black",getPalette(7)))+
    scale_size_discrete(range=c(4,6))+
    scale_alpha_discrete(range=c(.5,1))+ #(values=c(rep(0.2,7), 1))+
    guides(alpha=FALSE, size=FALSE) + #removes the legend 
    ggtitle(label = paste(title))+
    labs(x="", y=paste(units))+
    scale_x_discrete(labels = var_names) + # rev(unique(sort(dff$var))))+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=.9))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.title = element_blank())+
  
if(legend==FALSE){
  guides(color=FALSE) }
#pdf(file.path(figpath, "Fig2.pdf"), width = 7, height = 8)
  modfig <<- fig
}


### Loading models: 
load("/home/harold/Dropbox/gitfiles/germ/mod_time_pois.rdata") #this is my rstanarm  stanfit object.
load("/home/harold/Dropbox/gitfiles/germ/mod_rate.rdata") #this nd the next models are  rstanarm stanfit objects
load("/home/harold/Dropbox/gitfiles/germ/mod_gr.rdata")

# running function: 
fig1<-modplot(model=mod_time_pois, "b) Germination timing", "log(days)", legend=FALSE)
fig1<-  fig1+ theme(axis.text.x = element_blank())

fig2<-modplot(model=mod_rate, "a) Germination rate", "logit(fraction)", legend=TRUE) + theme(axis.text.x = element_blank())

fig3 <- modplot(model=mod_gr, "c) Growth rate", "cm/day", legend=FALSE)

## Putting them all together: 
fig2n<-fig2 + theme(legend.position ="none") 


grid<-plot_grid(fig2n, fig1, fig3, nrow=3, rel_heights = c(1,1,1.4))
#extract legend
legend <- get_legend( 
  # create some space to the left of the legend
#  fig2+ guides(color = guide_legend(nrow = 1)) +
 #   theme(legend.position = "bottom")
  fig2 + theme(legend.box.margin = margin(12, 0, 0, 12))
)
final  <- plot_grid(grid, legend, ncol=2, rel_widths = c(3,.8), axis='t', align="h")

ggsave("germ_figs_onepage.svg",final, device = "svg",width = 10,height=11, units="in")
#plot_grid(fig2n, fig1, fig3, ncol=1, align = "hv", axis="l",greedy=TRUE)
#plot_grid(fig2, fig1, fig3, ncol=1, rel_heights = c(1,.1,.7), align="v", axis = "l")

#saving the files: 

#png(file.path("germ_figs.png"), width = 7, height = 24, units="in", res=700)
#multiplot(fig1, fig2, fig3, cols=3)
#multiplot(fig2, fig1+theme(legend.position="none", axis.text.y=element_blank()))
#dev.off()
#pdf("germ_figs.pdf")
#fig1
#fig2
#fig3
#dev.off()
#save(fig1, fig2, fig3, file = "thesis_plots_raneff.Rdata")
