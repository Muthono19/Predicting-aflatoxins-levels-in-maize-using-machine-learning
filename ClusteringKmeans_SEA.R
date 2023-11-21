# ## initialization
# a script developed by Dr. Francis and modified by S. Gachoki
# create unique clusters to be used as unique space ID for
### Spatial cross-validation ###

### Load required libraries
library(RStoolbox)
library(raster)
library(rgdal)
library(sf)
library(lattice)
library(terra)
library(rasterVis)
library(RColorBrewer)

##Set the working directory
setwd("~/South Eastern Africa Aflatoxins")

### Load the precipitation long-term stacked layer
prec <- stack("clusterPrec_SEA.tif")
names(prec) <- c("Year2009",  "Year2010" , "Year2011" , "Year2012"  ,"Year2013" , "Year2014" , "Year2015" 
                ,"Year2016" , "Year2017"  ,"Year2018"  ,"Year2019" ,"Year2020" ,"Year2021", "Year2022")
## Calculate the average conditions
avgPrec <- mean(prec)
plot(avgPrec)

## Calculate how much each year deviates from the average conditions
year2009 <- prec[[1]] - avgPrec
year2010 <- prec[[2]] - avgPrec
year2011 <- prec[[3]] - avgPrec
year2012 <- prec[[4]] - avgPrec
year2013 <- prec[[5]] - avgPrec
year2014 <- prec[[6]] - avgPrec
year2015 <- prec[[7]] - avgPrec
year2016 <- prec[[8]] - avgPrec
year2017 <- prec[[9]] - avgPrec
year2018 <- prec[[10]] - avgPrec
year2019 <- prec[[11]] - avgPrec
year2020 <- prec[[12]] - avgPrec
year2021 <- prec[[13]] - avgPrec
year2022 <- prec[[14]] - avgPrec

deviation_average <- stack(year2009,  year2010 , year2011 , year2012  ,year2013 , year2014 , year2015 
                           ,year2016 , year2017  ,year2018  ,year2019 ,year2020 ,year2021, year2022)

names(deviation_average) <- c("Year2009",  "Year2010" , "Year2011" , "Year2012"  ,"Year2013" , "Year2014" , "Year2015" 
                 ,"Year2016" , "Year2017"  ,"Year2018"  ,"Year2019" ,"Year2020" ,"Year2021", "Year2022")

clorpalet <- colorRampPalette(brewer.pal(11,"RdYlBu"))
breaks <- c(-Inf,-500,-100, 0, 100, 200, 300, 400, 500,Inf)

png(file = "LongtermRainfallDeviations_SEA.png", width = 7000, height =7000, units = "px", res = 650, type = "cairo")
levelplot(deviation_average, at=breaks, col.regions=clorpalet(10))
dev.off()

## Extract values from the raster stack , calculate k-means clusters and write in one layer
kmeans.values <- getValues(prec)
kmeans.na <- which(!is.na(kmeans.values))
kmeans.values <- na.omit(kmeans.values)
kmeans.values.clust <- kmeans(kmeans.values,5,iter.max = 100,nstart = 2) # k=10
prec.2009 <- prec[[1]]
prec.2009[kmeans.na] <- kmeans.values.clust$cluster
plot(prec.2009)

writeRaster(prec.2009,"NewClusters_SEA.tif", overwrite=TRUE)

# ###**************************THE END********************####
