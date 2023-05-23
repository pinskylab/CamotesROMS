#2017 drifter visualization
library(maps)
library(mapdata)
library(maptools)
library(rgdal)
library(ggmap)
library(ggplot2)

drifters <- read.csv("empirical_data/2017drifters.csv")
map_new <- readOGR("leyte", "coastlines_leyte")
plot(map_new, xlim=c(124, 125 ), ylim=c(10.5,10.9))
par(new=TRUE)
plot(drifters$Latitude, drifters$Longitude, cex=0.1 )

ggplot() + geom_polygon(data = map_new, aes(x=long, y = lat, group = group), fill = "grey") + 
 coord_map(xlim=c(124, 125 ), ylim=c(10.5,11)) + geom_point(data = drifters, aes(x = Longitude, y = Latitude), color = "black", size = 1) 