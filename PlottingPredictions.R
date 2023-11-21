## Load the relevant libraries
library(raster)
library(rasterVis)
library(sf)
library(rgdal)
library(RColorBrewer)
library(lattice)
library(plotly)
library(maptools)
library(sp)
install.packages("maptools", repos="http://R-Forge.R-project.org")

### set the working directory
setwd("C:/Users/gstel/OneDrive/Desktop/IITA/Final analysis/Combined Final Analysis")

### Load the national boundaries
aoi <- readOGR("Kenya_Uganda_TanzaniaPRJ.shp")

## Load the predicted rasters
afla.pred <- stack("prediction_Farm_temp.tif","prediction_Future_Temp.tif","prediction_Farm_gap.tif","prediction_admin_Temp.tif")
names(afla.pred) <- c("","","","")

allnames <- c("D1","D1-2070","D2","D3")

myInverseColors <- rev(brewer.pal(11, "RdYlGn"))

myTheme=rasterTheme(region=myInverseColors)

png(file = "predictionsAll.png", width = 3500, height = 3500, units = "px", res = 700, type = "cairo")
levelplot(afla.pred,par.settings=myTheme, names.attr=allnames, layout=c(2,2), margin=FALSE,ylab="",
          panel=function(...) {
            panel.levelplot(...)
            sp.polygons(aoi, fill = NA, col = "black", lwd=0.3)
          })
dev.off()


?levelplot

