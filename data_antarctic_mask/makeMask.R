# Make a mask incorporating land as well as Antarctic ice sheets
# Ryan R Reisinger
# 2018-04-18

setwd("~/RAATD_01/RAATD/antarcticMask")

library(raster)
library(raadtools)

# Load the coast mask
# Obtained from:
# https://nsidc.org/data/nsidc-0709/versions/2#
# on 2018/04/17

coast <- raster("Mask_Antarctica_v02.tif")
coast[coast < 50] <- NA

# Make an exmpty raster
lims <- c(-180, 180, -80, -40)
grd <- raster(extent(lims), resolution = c(0.1, 0.1), crs = "+proj=longlat +datum=WGS84")

# Project the mask
coastProj <- projectRaster(from = coast, to = grd)

# Fix the land blip
rst <- coastProj
temp <- matrix(values(rst), nrow = dim(rst)[1], byrow = TRUE)
temp[380:400, 1:13] <- minValue(rst)
values(rst) <- temp
coastProj <- rst

# Fill in NA values
coastProj[is.na(coastProj)] <- 0

# Create the mask
coastProj[coastProj > 0] <- NA

# Additional land mask
bath <- readderivaadc("bathymetry", xylim = lims)
bath <- resample(bath, grd, method = "bilinear")
bath <- crop(bath, extent(grd))
bath[bath >= 0] <- NA

# Combine the masks
msk <- mask(coastProj, bath)

# Write output
writeRaster(msk, "mask.grd", format = "raster", overwrite = T)
