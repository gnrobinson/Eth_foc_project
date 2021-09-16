########################################################################################
# Title: Plotting isolate collection across Ethiopia
# Author: Guy Robinson
########################################################################################
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
# https://www.r-bloggers.com/2018/08/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/

#loading required libraries
library(tidyverse)
library(sf)
library(ggplot2)
library(ggspatial)
library(rgeos)
library(rnaturalearth)
library(rgdal)
library(raster)

#load fungal isolate collection data
setwd("/Path/to/Directory/")
geog <- read_delim(file = "geography_data.txt", col_names = TRUE, delim = "\t")
geog$mat <- as.character(geog$mat)

##### Data structure for analysis ####
#isolates   latitude longitude   mat `sub-groups``
# <chr>         <dbl>     <dbl> <chr> <chr>       
# 1 Fo-Et-0000     9.42      42.0     1 T           
# 2 Fo-Et-0001     9.42      41.6     2 Sing        
# 3 Fo-Et-0002     8.95      39.1     2 Sing        
# 4 Fo-Et-0003     9.01      39.5     2 Sing        
# 5 Fo-Et-0008     9.98      39.2     2 Sing        
# 6 Fo-Et-0011     7.81      39.1     1 Sing  



#Download world map data from rnaturalearth
world <- ne_countries(scale = "medium", returnclass = "sf")

#Download river data via rgdal-- Dataset downloaded from naturalearthdata.com
rivers <- readOGR('ne_50m_rivers_lake_centerlines.shp',
                     'ne_50m_rivers_lake_centerlines')

#Download altitude data for plotting mountain ranges
dem.raster <- getData("alt", country=c('ETH'), download = TRUE)
dem.m  <-  rasterToPoints(dem.raster) 
dem.df <-  data.frame(dem.m) 
colnames(dem.df) = c("lon", "lat", "alt")

slope.raster <- terrain(dem.raster, opt='slope') #
aspect.raster <- terrain(dem.raster, opt='aspect') #
hill.raster <- hillShade(slope.raster, aspect.raster, 40, 270) #
hill.m <- rasterToPoints(hill.raster) #
hill.df <-  data.frame(hill.m) #
colnames(hill.df) <- c("lon", "lat", "hill")
hill.df$log_hill <- abs(log10(hill.df$hill))
hill.df$inv_log <- abs(1/log10(hill.df$hill))

#adding colors for plotting
cols = c(
  "#b8dedb",
  "#b8dedb",
  "#b8dedb",
  "#b8dedb",
  "#b8dedb",
  "#9cb600",
  "#b8dedb",
  "#9cb600",
  "#37d279",
  "#b8dedb",
  "#ff40ff",
  "#ff40ff",
  "#37d279",
  "#bf8a93",
  "#074FFF",
  "#b8dedb",
  "#074FFF",
  "#3ab89a",
  "#3ab89a",
  "#074FFF",
  "#b8dedb",
  "#ff2600",
  "#b8dedb",
  "#6ecbe3",
  "#074FFF",
  "#9cb600",
  "#ff40ff",
  "#b8dedb",
  "#ff2600",
  "#9cb600",
  "#074FFF",
  "#074FFF",
  "#ff2600",
  "#bf8a93",
  "#bf8a93",
  "#b8dedb",
  "#b8dedb",
  "#bb9d43",
  "#bb9d43",
  "#b8dedb",
  "#37d279",
  "#b8dedb",
  "#b8dedb",
  "#9cb600",
  "#ff2600",
  "#9cb600",
  "#ff2600",
  "#d84800",
  "#b8dedb",
  "#b8dedb",
  "#b8dedb",
  "#37d279",
  "#d84800",
  "#b8dedb",
  "#b8dedb",
  "#9e9873",
  "#b8dedb",
  "#b8dedb",
  "#007d7d",
  "#b8dedb",
  "#b8dedb",
  "#ff40ff",
  "#b8dedb",
  "#d84800",
  "#37d279",
  "#ff40ff",
  "#d84800",
  "#b8dedb",
  "#b8dedb",
  "#d84800",
  "#d84800",
  "#d84800",
  "#b8dedb",
  "#ff40ff",
  "#b8dedb",
  "#b8dedb",
  "#d84800",
  "#37d279",
  "#942092",
  "#3ab89a",
  "#d84800",
  "#007d7d",
  "#fffc00",
  "#b8dedb",
  "#377eb8",
  "#fffc00",
  "#ffcee4",
  "#b8dedb",
  "#fffc00",
  "#b8dedb",
  "#7f9fe6",
  "#7f9fe6",
  "#37d279",
  "#fffc00",
  "#fffc00",
  "#ff2600",
  "#d84800",
  "#b8dedb",
  "#37d279",
  "#9cb600",
  "#6ecbe3",
  "#007d7d",
  "#37d279",
  "#b8dedb",
  "#37d279",
  "#b8dedb",
  "#ffcee4",
  "#377eb8",
  "#377eb8",
  "#3ab89a",
  "#37d279",
  "#ff40ff",
  "#ff2600",
  "#007d7d",
  "#9cb600",
  "#6ecbe3",
  "#d84800",
  "#d84800",
  "#b8dedb",
  "#37d279",
  "#942092",
  "#9e9873",
  "#007d7d",
  "#b8dedb",
  "#ff2600",
  "#ff2600",
  "#d84800",
  "#b8dedb",
  "#074FFF",
  "#d84800",
  "#37d279",
  "#37d279",
  "#b8dedb")

colnames(geog)[4] <- "Mat Locus"

#plotting geographic data
ggplot(data = world) +
  geom_sf(fill = 'white', lwd = 0.1) +
  coord_sf(xlim = c(32.00, 47.00), ylim = c(2.00, 18.50), expand = FALSE) +
  geom_raster(data = hill.df, aes(lon, lat, fill = hill)) +
  geom_sf(fill = NA, lwd = 0.1) +
  coord_sf(xlim = c(32.00, 47.00), ylim = c(2.00, 18.50), expand = FALSE) +
  geom_point(data = geog, aes(x=longitude, y=latitude, shape = `Mat Locus`), size = 5, color = cols) +
  scale_fill_gradient2(low = '#ebf2f2', high = '#0d0f0f', midpoint = 0.63) +
  scale_shape_manual(values=c(3, 1))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = 'grey92'),
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21)) + #panel.grid.major = element_line(color = 'grey', linetype = 'dotted', size = 0.5)
  annotate(geom = 'text', x = 42, y = 7, label = 'Ethiopia', fontface = 'italic', color = 'grey22', size = 10) +
  annotate(geom = 'text', x = 38, y = 16, label = 'Eritrea', fontface = 'italic', color = 'grey22', size = 10) +
  annotate(geom = 'text', x = 42.5, y = 12, label = 'Djibouti', fontface = 'italic', color = 'grey22', size = 10) +
  annotation_scale(location = 'bl', width_hint = 0.5, height = unit(0.4, "cm"), text_cex = 2) +
  geom_path(data=rivers, aes(x = long, y = lat, group = group), color = '#a2baeb')
  
