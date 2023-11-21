library(raster)
library(rgdal)


setwd("~/South Eastern Africa Aflatoxins")

aoi <- readOGR("AOI_Boundaries.shp")
extent <- extent(aoi)
resolution <- 0.01  # 1 km resolution
plot(aoi)

# Calculate the number of cells in the longitude raster
ncells_longitude <- ceiling((extent@xmax - extent@xmin) / resolution)

# Create a Raster for Longitudes
longitude_raster <- raster(extent, ncol = ncells_longitude)

# Generate Longitudes
longitude_values <- seq(extent@xmin, extent@xmax, by = resolution)

# Set Raster Values for Longitudes
longitude_raster[] <- longitude_values

plot(longitude_raster)

# Write the Longitude Raster to a GeoTIFF File
writeRaster(longitude_raster, "Longitudes.tif", format = "GTiff", overwrite=TRUE)
