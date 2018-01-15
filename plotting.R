## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(shinystan.rstudio = TRUE)
options(mc.cores = parallel::detectCores())

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/undergrads/harold/analyses/germination_stan") 
} else 
  setwd("C:/Users/Owner/Documents/GitHub/germination_stan")

source("http://peterhaschke.com/Code/multiplot.R") #so that the multiplot function works 
figpath<-"C:/Users/Owner/Documents/Github/germination_stan"


#vertical position dodge 
# see https://raw.githubusercontent.com/jaredlander/coefplot/master/R/position.r
#  https://www.jaredlander.com/2013/02/vertical-dodging-in-ggplot2/ 
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/jaredlander/coefplot/master/R/position.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

##libraries
library(rstan)
library(shinystan)
library(lme4)
library(ggplot2)
library(rstanarm)
library(brms)
library(tidyr)

####### Days from Germination Plot:
#######
#######



#data
load("C:/Users/Owner/Documents/Thesis/Stan/mod_time_pois_brm.Rdata")

m<-mod_time_pois_brm
m.int<-posterior_interval(m)
sum.m<-summary(m)
cri.f<-as.data.frame(sum.m$fixed[,c("Estimate", "l-95% CI", "u-95% CI")])
cri.f<-cri.f[-1,] #removing the intercept 
fdf1<-as.data.frame(rbind(as.vector(cri.f[,1]), as.vector(cri.f[,2]), as.vector(cri.f[,3])))
fdf2<-cbind(fdf1, c(0, 0, 0) , c("Estimate", "2.5%", "95%"))
names(fdf2)<-c(rownames(cri.f), "sp", "perc")

cri.r<-(ranef(m, summary = TRUE, robust = FALSE,
     probs = c(0.025, 0.975)))$sp
cri.r2<-cri.r[, ,-1]
cri.r2<-cri.r2[,-2,]
dims<-dim(cri.r2)
twoDimMat <- matrix(cri.r2, prod(dims[1:2]), dims[3])
mat2<-cbind(twoDimMat, c(rep(1:7, length.out=21)), rep(c("Estimate", "2.5%", "95%"), each=7))
df<-as.data.frame(mat2)
names(df)<-c(rownames(cri.f), "sp", "perc")
dftot<-rbind(fdf2, df)
dflong<- gather(dftot, var, value, origin:`origin:strat:temp3`, factor_key=TRUE)

#adding the coef estiamtes to the random effect values 
for (i in seq(from=1,to=nrow(dflong), by=24)) {
  for (j in seq(from=3, to=23, by=1)) {
    dflong$value[i+j]<- as.numeric(dflong$value[i+j]) + as.numeric(dflong$value[i])
  }
}
dflong$rndm<-ifelse(dflong$sp>0, 2, 1)
dfwide<-spread(dflong, perc, value)
dfwide[,4:6] <- as.data.frame(lapply(c(dfwide[,4:6]), as.numeric ))
dfwide$sp<-as.factor(dfwide$sp)
## plotting

pdf(file.path(figpath, "Fig1.pdf"), width = 7, height = 8)
pd <- position_dodgev(height = -0.5)


fig1 <-ggplot(dfwide, aes(x=Estimate, y=var, color=factor(sp), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`95%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels = c("Fixed effects", "CAPBUR", "CHEMAJ", "DACGLO", "PLALAN", "PLAMAJ", "RUMCRI", "TAROFF"),
                      values=c("blue", "red", "orangered1", "sienna4", "green4", "green1", "purple2", "magenta2"))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Days to Germination")+ 
  scale_y_discrete(limits = rev(unique(sort(dfwide$var))))
fig1
dev.off()

#####
#####
#####Germination Rate Plot 

#data: 
load("C:/Users/Owner/Documents/Thesis/Stan/mod_rate.Rdata")

m<-mod_rate
sum.m<-summary(m, pars=c("(Intercept)", "origin", "strat", "temp1", "temp2", "temp3", "origin:strat", "origin:temp1", "origin:temp2",
 "origin:temp3", "strat:temp1", "strat:temp2", "strat:temp3", "origin:strat:temp1", "origin:strat:temp2", "origin:strat:temp3"
, "b[origin sp:1]", "b[strat sp:1]", "b[temp1 sp:1]", "b[temp2 sp:1]", "b[temp3 sp:1]", "b[origin:strat sp:1]",
"b[origin:temp1 sp:1]", "b[origin:temp2 sp:1]", "b[origin:temp3 sp:1]", "b[strat:temp1 sp:1]", "b[strat:temp2 sp:1]", "b[strat:temp3 sp:1]", "b[origin:strat:temp1 sp:1]",
  "b[origin:strat:temp2 sp:1]", "b[origin:strat:temp3 sp:1]"
,"b[origin sp:2]", "b[strat sp:2]", "b[temp1 sp:2]", "b[temp2 sp:2]", "b[temp3 sp:2]", "b[origin:strat sp:2]",
                                            "b[origin:temp1 sp:2]", "b[origin:temp2 sp:2]", "b[origin:temp3 sp:2]","b[strat:temp1 sp:2]", "b[strat:temp2 sp:2]", "b[strat:temp3 sp:2]", "b[origin:strat:temp1 sp:2]", 
                                            "b[origin:strat:temp2 sp:2]", "b[origin:strat:temp3 sp:2]"
,"b[origin sp:3]", "b[strat sp:3]", "b[temp1 sp:3]", "b[temp2 sp:3]", "b[temp3 sp:3]", "b[origin:strat sp:3]",
                                            "b[origin:temp1 sp:3]", "b[origin:temp2 sp:3]",  "b[origin:temp3 sp:3]", "b[strat:temp1 sp:3]", "b[strat:temp2 sp:3]", "b[strat:temp3 sp:3]", "b[origin:strat:temp1 sp:3]",
                                            "b[origin:strat:temp2 sp:3]", "b[origin:strat:temp3 sp:3]"
,"b[origin sp:4]", "b[strat sp:4]", "b[temp1 sp:4]", "b[temp2 sp:4]", "b[temp3 sp:4]", "b[origin:strat sp:4]",
"b[origin:temp1 sp:4]", "b[origin:temp2 sp:4]",  "b[origin:temp3 sp:4]", "b[strat:temp1 sp:4]", "b[strat:temp2 sp:4]", "b[strat:temp3 sp:4]", "b[origin:strat:temp1 sp:4]",
 "b[origin:strat:temp2 sp:4]", "b[origin:strat:temp3 sp:4]"
,"b[origin sp:5]", "b[strat sp:5]", "b[temp1 sp:5]", "b[temp2 sp:5]", "b[temp3 sp:5]", "b[origin:strat sp:5]",
"b[origin:temp1 sp:5]", "b[origin:temp2 sp:5]",  "b[origin:temp3 sp:5]", "b[strat:temp1 sp:5]", "b[strat:temp2 sp:5]", "b[strat:temp3 sp:5]", "b[origin:strat:temp1 sp:5]",
 "b[origin:strat:temp2 sp:5]", "b[origin:strat:temp3 sp:5]"
,"b[origin sp:6]", "b[strat sp:6]", "b[temp1 sp:6]", "b[temp2 sp:6]", "b[temp3 sp:6]", "b[origin:strat sp:6]",
"b[origin:temp1 sp:6]", "b[origin:temp2 sp:6]",  "b[origin:temp3 sp:6]","b[strat:temp1 sp:6]", "b[strat:temp2 sp:6]", "b[strat:temp3 sp:6]", "b[origin:strat:temp1 sp:6]",
"b[origin:strat:temp2 sp:6]", "b[origin:strat:temp3 sp:6]"
,"b[origin sp:7]", "b[strat sp:7]", "b[temp1 sp:7]", "b[temp2 sp:7]", "b[temp3 sp:7]", "b[origin:strat sp:7]",
 "b[origin:temp1 sp:7]", "b[origin:temp2 sp:7]",  "b[origin:temp3 sp:7]","b[strat:temp1 sp:7]", "b[strat:temp2 sp:7]", "b[strat:temp3 sp:7]", "b[origin:strat:temp1 sp:7]",
"b[origin:strat:temp2 sp:7]", "b[origin:strat:temp3 sp:7]"))

cri.f<-sum.m[c(2:121),c(1,4,8)] #just selecting the mean and 95% CI. Removing the intercept 
fdf<-data.frame(cri.f)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rev(rownames(fdf)[1:15]), each=7))),
                             as.numeric(as.character(fdf$mean)), 
                             as.numeric(as.character(fdf$X2.5.)), 
                             as.numeric(as.character(fdf$X97.5.)), 
                             as.numeric(c(rep(1, 15), rep(2, 105))), 
                            as.numeric( c(rep(0,15), rep(seq(1:7),15 )))))

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "sp")

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Fixed effect estimates 
fixed<-c(rep(0, 15), rep(as.numeric(rev(fdf2[c(1:15),2])), each=7))
dff<-fdf2
#adding the fixef estiamtes to the random effect values 

dff$Estimate<-fdf2$Estimate+fixed
dff$`2.5%`<-fdf2$`2.5%`+fixed
dff$`97.5%`<-fdf2$`97.5%`+fixed

## plotting

pdf(file.path(figpath, "Fig2.pdf"), width = 7, height = 8)
pd <- position_dodgev(height = -0.5)


fig2<-ggplot(dff, aes(x=Estimate, y=var, color=factor(sp), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels = c("Fixed effects", "CAPBUR", "CHEMAJ", "DACGLO", "PLALAN", "PLAMAJ", "RUMCRI", "TAROFF"),
                      values=c("blue", "red", "orangered1", "sienna4", "green4", "green1", "purple2", "magenta2"))+  scale_alpha_manual(values=c(1, 0.5))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Germination rate")+
  scale_y_discrete(limits = rev(unique(sort(dff$var))))
fig2
dev.off()

#####
#####
#####
#next, growth rate:

#data: 
load("C:/Users/Owner/Documents/Thesis/Stan/mod_gr.Rdata")

m<-mod_gr
sum.m<-summary(m, pars=c("(Intercept)", "origin", "strat", "temp1", "temp2", "temp3", "origin:strat", "origin:temp1", "origin:temp2",
                         "origin:temp3", "strat:temp1", "strat:temp2", "strat:temp3", "origin:strat:temp1", "origin:strat:temp2", "origin:strat:temp3"
                         , "b[origin sp:1]", "b[strat sp:1]", "b[temp1 sp:1]", "b[temp2 sp:1]", "b[temp3 sp:1]", "b[origin:strat sp:1]",
                         "b[origin:temp1 sp:1]", "b[origin:temp2 sp:1]", "b[origin:temp3 sp:1]", "b[strat:temp1 sp:1]", "b[strat:temp2 sp:1]", "b[strat:temp3 sp:1]", "b[origin:strat:temp1 sp:1]",
                         "b[origin:strat:temp2 sp:1]", "b[origin:strat:temp3 sp:1]"
                         ,"b[origin sp:2]", "b[strat sp:2]", "b[temp1 sp:2]", "b[temp2 sp:2]", "b[temp3 sp:2]", "b[origin:strat sp:2]",
                         "b[origin:temp1 sp:2]", "b[origin:temp2 sp:2]", "b[origin:temp3 sp:2]","b[strat:temp1 sp:2]", "b[strat:temp2 sp:2]", "b[strat:temp3 sp:2]", "b[origin:strat:temp1 sp:2]", 
                         "b[origin:strat:temp2 sp:2]", "b[origin:strat:temp3 sp:2]"
                         ,"b[origin sp:3]", "b[strat sp:3]", "b[temp1 sp:3]", "b[temp2 sp:3]", "b[temp3 sp:3]", "b[origin:strat sp:3]",
                         "b[origin:temp1 sp:3]", "b[origin:temp2 sp:3]",  "b[origin:temp3 sp:3]", "b[strat:temp1 sp:3]", "b[strat:temp2 sp:3]", "b[strat:temp3 sp:3]", "b[origin:strat:temp1 sp:3]",
                         "b[origin:strat:temp2 sp:3]", "b[origin:strat:temp3 sp:3]"
                         ,"b[origin sp:4]", "b[strat sp:4]", "b[temp1 sp:4]", "b[temp2 sp:4]", "b[temp3 sp:4]", "b[origin:strat sp:4]",
                         "b[origin:temp1 sp:4]", "b[origin:temp2 sp:4]",  "b[origin:temp3 sp:4]", "b[strat:temp1 sp:4]", "b[strat:temp2 sp:4]", "b[strat:temp3 sp:4]", "b[origin:strat:temp1 sp:4]",
                         "b[origin:strat:temp2 sp:4]", "b[origin:strat:temp3 sp:4]"
                         ,"b[origin sp:5]", "b[strat sp:5]", "b[temp1 sp:5]", "b[temp2 sp:5]", "b[temp3 sp:5]", "b[origin:strat sp:5]",
                         "b[origin:temp1 sp:5]", "b[origin:temp2 sp:5]",  "b[origin:temp3 sp:5]", "b[strat:temp1 sp:5]", "b[strat:temp2 sp:5]", "b[strat:temp3 sp:5]", "b[origin:strat:temp1 sp:5]",
                         "b[origin:strat:temp2 sp:5]", "b[origin:strat:temp3 sp:5]"
                         ,"b[origin sp:6]", "b[strat sp:6]", "b[temp1 sp:6]", "b[temp2 sp:6]", "b[temp3 sp:6]", "b[origin:strat sp:6]",
                         "b[origin:temp1 sp:6]", "b[origin:temp2 sp:6]",  "b[origin:temp3 sp:6]","b[strat:temp1 sp:6]", "b[strat:temp2 sp:6]", "b[strat:temp3 sp:6]", "b[origin:strat:temp1 sp:6]",
                         "b[origin:strat:temp2 sp:6]", "b[origin:strat:temp3 sp:6]"
                         ,"b[origin sp:7]", "b[strat sp:7]", "b[temp1 sp:7]", "b[temp2 sp:7]", "b[temp3 sp:7]", "b[origin:strat sp:7]",
                         "b[origin:temp1 sp:7]", "b[origin:temp2 sp:7]",  "b[origin:temp3 sp:7]","b[strat:temp1 sp:7]", "b[strat:temp2 sp:7]", "b[strat:temp3 sp:7]", "b[origin:strat:temp1 sp:7]",
                         "b[origin:strat:temp2 sp:7]", "b[origin:strat:temp3 sp:7]"))

cri.f<-sum.m[c(2:121),c(1,4,8)] #just selecting the mean and 95% CI. Removing the intercept 
fdf<-data.frame(cri.f)
#binding 
fdf2<-as.data.frame(
  cbind(
    (c(rownames(fdf)[1:15], rep(rev(rownames(fdf)[1:15]), each=7))),
    as.numeric(as.character(fdf$mean)), 
    as.numeric(as.character(fdf$X2.5.)), 
    as.numeric(as.character(fdf$X97.5.)), 
    as.numeric(c(rep(1, 15), rep(2, 105))), 
    as.numeric( c(rep(0,15), rep(seq(1:7),15 )))))

names(fdf2)<-c("var", "Estimate", colnames(cri.f)[c(2,3)], "rndm", "sp")

fdf2$Estimate<-as.numeric(fdf2$Estimate)      
fdf2$`2.5%`<-as.numeric(fdf2$`2.5%`)
fdf2$`97.5%`<-as.numeric(fdf2$`97.5%`)      


#Fixed effect estimates 
fixed<-c(rep(0, 15), rep(as.numeric(rev(fdf2[c(1:15),2])), each=7))
dff<-fdf2
#adding the fixef estiamtes to the random effect values 

dff$Estimate<-fdf2$Estimate+fixed
dff$`2.5%`<-fdf2$`2.5%`+fixed
dff$`97.5%`<-fdf2$`97.5%`+fixed

## plotting

pdf(file.path(figpath, "Fig3.pdf"), width = 7, height = 8)
pd <- position_dodgev(height = -0.5)


fig3<-ggplot(dff, aes(x=Estimate, y=var, color=factor(sp), size=factor(rndm), alpha=factor(rndm)))+
  geom_point(position =pd, size=4)+
  geom_errorbarh(aes(xmin=(`2.5%`), xmax=(`97.5%`)), position=pd, size=.5, height =0)+
  geom_vline(xintercept=0)+
  scale_colour_manual(labels = c("Fixed effects", "CAPBUR", "CHEMAJ", "DACGLO", "PLALAN", "PLAMAJ", "RUMCRI", "TAROFF"),
                      values=c("blue", "red", "orangered1", "sienna4", "green4", "green1", "purple2", "magenta2"))+  scale_alpha_manual(values=c(1, 0.5))+
  scale_shape_manual(labels="", values=c("1"=16,"2"=16))+
  scale_alpha_manual(values=c(1, 0.5))+
  guides(alpha=FALSE) + #removes the legend 
  ggtitle(label = "Growth rate")+
  scale_y_discrete(limits = rev(unique(sort(dff$var))))
fig3
dev.off()
png(file.path(figpath, "Eyster_figs.png"), width = 7, height = 24, units="in", res=200)
multiplot(fig1, fig2, fig3)
dev.off()






# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    
    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    
    #   # Suppress warning message since it's not reliable
    #     if (!zero_range(range(heights))) {
    #       warning(name, " requires constant height: output may be incorrect",
    #         call. = FALSE)
    #     }
    height <- heights[1]
  }
  
  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]
  
  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }
  
  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
  if (nrow(df) == 1) return(df)
  
  n <- nrow(df) + 1
  x <- ifelse(is.na(df$x), 0, df$x)
  if (all(is.na(df$y))) {
    heights <- rep(NA, n)
  } else {
    heights <- c(0, cumsum(x))
  }
  
  df$xmin <- heights[-n]
  df$xmax <- heights[-1]
  df$x <- df$xmax
  df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
  stacked <- pos_stackv(df, height)
  stacked$xmin <- stacked$xmin / max(stacked$xmax)
  stacked$xmax <- stacked$xmax / max(stacked$xmax)
  stacked$x <- stacked$xmax
  stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)
  
  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }
  
  d_height <- max(df$ymax - df$ymin)
  
  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()
  
  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group)))
  
  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  
  df
}


#' Adjust position by dodging overlaps to the side.
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
#' @examples
#' ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#'   geom_bar(position = "dodge")
#'
#' ggplot(diamonds, aes(price, fill = cut)) +
#'   geom_histogram(position="dodge")
#' # see ?geom_boxplot and ?geom_bar for more examples
#'
#' # To dodge items with different heights, you need to be explicit
#' df <- data.frame(x=c("a","a","b","b"), y=2:5, g = rep(1:2, 2))
#' p <- ggplot(df, aes(x, y, group = g)) +
#'   geom_bar(
#'     stat = "identity", position = "dodge",
#'     fill = "grey50", colour = "black"
#'   )
#' p
#'
#' # A line range has no height:
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1), position = "dodge")
#' # You need to explicitly specify the height for dodging
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1),
#'   position = position_dodge(width = 0.9))
#'
#' # Similarly with error bars:
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1), width = 0.2,
#'   position = "dodge")
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1, height = 0.2),
#'   position = position_dodge(width = 0.90))
#'
position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}



PositionDodgeV <- ggproto("PositionDodgeV", Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                            if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                              warning("height not defined. Set with `position_dodgev(height = ?)`",
                                      call. = FALSE)
                            }
                            list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                            collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)

