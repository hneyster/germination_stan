library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tools)
library(maps)
library(rgeos)
library(cowplot)
library(RColorBrewer)
#ref: https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-2.html

getPalette = colorRampPalette(brewer.pal(13, "Set1"))


clim_raw<-read.csv("/home/harold/github/germination_stan/avgtemps.csv")
load("/home/harold/github/germination_stan/germs.Rdata") #cleaned and processed real data
germs.y<-(subset(germs, germinated==1 & sp!="PLAMED" &  sp!="PLACOR"))    #just the data from seeds that germianted, and taking out the congenerics
loc<-as.numeric(as.factor(germs.y$loc))
locs_names <- cbind(germs.y$loc, loc) %>% as.data.frame() %>% unique() 
locs_names<-locs_names[order(as.numeric(as.character(locs_names$loc))),]

clim<-merge(clim_raw, locs_names, by.x="loc_name", by.y="V1", all = FALSE) 
names(clim)[7] <- "loc_num"

clim_us<-clim[clim$lon<0,]
clim_eu <- clim[clim$lon>0,]
#jittering the two very-close point in Austria 

sites_eu <- st_as_sf(clim_eu, coords = c("lon", "lat"), 
                   crs = 4326, agr = "constant")
sites_us <-st_as_sf(clim_us, coords = c("lon", "lat"), 
                    crs = 4326, agr = "constant")

sites_eu[clim_eu$loc_name=="Europe--Zell Am See, Austria","lon"]<-  
  clim_eu[clim_eu$loc_name=="Europe--Zell Am See, Austria","lon"] + .2

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(as.character(states$ID))

us<-ggplot(data = world) + geom_sf() + geom_sf(data=sites_us, size=5, color=getPalette(13)[11:13], shape=17)  +
  geom_sf(data = states, fill = NA) +
coord_sf(xlim = c(-75, -70), ylim = c(40, 45))+
  theme(panel.grid.major = element_line(colour = "white")) + theme_bw()
#pdf("us_collections.pdf", width = 6)
 # geom_text(data = states[states$ID=="Massachusetts",], aes(X, Y, label = ID), size = 4) 

eu<-ggplot(data = world) + geom_sf() + geom_sf(data=sites_eu, size=5, color=getPalette(16)[1:13])  +
  geom_sf(data = states, fill = NA) +
  coord_sf(xlim = c(-5, 20), ylim = c(40, 60))+
  theme(panel.grid.major = element_line(colour = "white")) + theme_bw()

#p1<-plot_grid(us, eu, align="h", ncol=2, labels='AUTO')
#ggsave(filename = "collections.pdf",plot = p1, width = 6, height = 6, units="in","pdf")


#Adding climate plot, now with lines: 
al<-.8
clim$cont <- factor(as.numeric(c(rep(0,13),rep(1,3))))

climwide <- clim
names(climwide)[names(climwide) %in% c('Mar',"Apr","May")] <- c('temp.1','temp.2','temp.3')
climlong <- reshape(climwide, direction="long",  varying=4:6)
climlong$cont <- ifelse(climlong$cont==0, 'Europe','North America')
clim_plot <- ggplot(data = climlong, aes(x = time, y = temp, color = factor(as.numeric(climlong$loc_num)),shape=climlong$cont), alpha=al) + 
  geom_line(position=position_dodge(-.4))+ geom_point(size = 4,position=position_dodge(-.4), alpha = .5)+ facet_grid(.~climlong$cont) + theme_bw() +   
  scale_color_manual(labels=clim$short, values=getPalette(nrow(clim)))+
  scale_x_continuous(breaks=c(1,2,3), labels=c("Mar","Apr","May")) + 
  ylab(expression(paste("Mean temperature 1970-2000 (",~degree*C,")")))+
  xlab("Month")+theme_bw() + theme(legend.position = "none")

#ggsave(filename = "climate_plot.pdf",plot = clim_plot, width = 6, height = 6, units="in","pdf")#Making it easier to differentiate continents: 
p0 <- plot_grid(eu,us)
p1<-plot_grid(clim_plot,p0, ncol = 1)
ggsave(filename = "sampling_sites.pdf",plot = p1, width = 6, height = 6, units="in","pdf")


#p1<-plot_grid(eu, us, clim_plot,align="hv", ncol=3, labels='auto',axis = 'tblr',rel_widths = c(1,1,2), rel_heights = c(1,1,1))
#clim_plot<-  ggplot(data=clim) + 
#geom_point(position=position_dodge(-.4), aes(y=clim$Mar, x=1, color=factor(as.numeric(clim$loc_num)),shape=clim$cont), size=5, alpha=al)+
#  geom_point(position=position_dodge(-.4),aes(x=2, y=clim$Apr, color=factor(as.numeric(clim$loc_num)),shape=clim$cont), size=5, alpha=al)+
#  geom_point(position=position_dodge(-.4),aes(x=3, y=clim$May, color=factor(as.numeric(clim$loc_num)),shape=clim$cont), size=5, alpha=al)+
#  scale_color_manual(labels=clim$short, values=getPalette(nrow(clim)))+
#  scale_x_continuous(breaks=c(1,2,3), labels=c("Mar","Apr","May")) + 
#  ylab(expression(paste("Mean temperature 1970-2000 (",~degree*C,")")))+
# xlab("Month")+theme_bw() + theme(legend.position = "none")
