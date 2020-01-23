## Started 6 Jan 2020 ##
## By Lizzie Wolkovich ##
## Modified 20 Jan 2020 by Harold Eyster ##
## Pulling climate data ##

# Help from ...
# Here: https://www.benjaminbell.co.uk/2018/01/extracting-data-and-making-climate-maps.html
# Cat and JD

# Housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)
setwd("~/Documents/GitHub/germination_stan")

# Load libraries
library(lubridate)
library(raster)

# Get the data
d_us <- read.csv("US_lat_long.csv")
d <- read.csv("europe_lat_long.csv")
d<-rbind(d,d_us)
dsm <- d[1:10,]
names(d) <- c("lat", "lon", "loc_name")
d <- d[c("lon", "lat", "loc_name")] # Stupidly you need them in this order to use extract...
# below from: http://worldclim.org/version2
temp3 <- raster("~/Documents/thesis/climate data/wc2.0_30s_tavg/wc2.0_30s_tavg_03.tif") # March
temp4 <- raster("~/Documents/thesis/climate data/wc2.0_30s_tavg/wc2.0_30s_tavg_04.tif") # Apr
temp5 <- raster("~/Documents/thesis/climate data/wc2.0_30s_tavg/wc2.0_30s_tavg_05.tif") # May

# Plot to check ... 
#plot(temp3[[1]], xlim=c(0,50), ylim=c(25, 75))
#points(lat~lon, data=d, cex=0.5, pch=16)

# Create a data.frame with sample site coordinates
temp.data <- d
temp.data$Mar <- extract(temp3, d[,1:2])
temp.data$Apr <- extract(temp4, d[,1:2])
temp.data$May <- extract(temp5, d[,1:2])

write.csv(temp.data, "avgtemps.csv", row.names=FALSE)
