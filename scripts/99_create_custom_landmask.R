## Create custom landmask

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

library(raster)

r <- raster("./data_landmask/ETOPO1_Ice_c_geotiff.tif")

these_which <- coordinates(r)[,2] > -40 # Which cells are outside the study area
these_values <- values(r) # Get all raster values
these_values[these_which] <- 1000 # Set those values outside study area as land
r <- setValues(r, these_values) # Reassign the raster values

# Write new landmask
writeRaster(r, "./data_landmask/ETOPO1_custom.tif")
