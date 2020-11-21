## Compile environmental variables

setwd("/mnt/home/ryan/humpbacks/megaPrediction")

library(raster)
library(raadtools)
library(grec)

#------------------------------
## Set up raster
lims <- c(-180, +180, -80, -40)
lims_pad <- c(-180, +180, -80, -30) # Where edge effects need to be avoided
wgs <- "+proj=longlat +datum=WGS84 +no_defs"
r <- raster(ext = extent(lims), resolution = c(0.1, 0.1), crs = wgs)
r_pad <- raster(ext = extent(lims_pad), resolution = c(0.1, 0.1), crs = wgs)

## The same raster on 0-360, for some products
r_360 <- raster(ext = extent(0, 360, -80, -40), resolution = c(0.1, 0.1), crs = wgs)
r_360_pad <- raster(ext = extent(0, 360, -80, -30), resolution = c(0.1, 0.1), crs = wgs)

#------------------------------
# Dates
# Only summer: November - March
which_months <- c("11", "12", "01", "02", "03")
dates <- data.frame(dates = seq(as.POSIXct("1999-11-01"), as.POSIXct("2019-03-31"), 86400*7))
dates$month <- format(dates$dates, "%m")
dates<- dates[dates$month %in% which_months, ]
# Austral year for calculating intra-season variance
dates$year <- format(dates$dates, "%Y")
dates$austral_year <- dates$year
dates[dates$month %in% c("01", "02", "03"), ]$austral_year <- as.integer(dates[dates$month %in% c("01", "02", "03"), ]$austral_year) - 1


#------------------------------
# Masks

## load mask already created (see folder 'antarcticMask')
## this mask includes the ice-shelves
landmask <- raster("./data_antarctic_mask/mask.grd")

## ice mask
ice <- mean(readice(dates$dates, setNA = T, latest = F), na.rm = T)
rd3dummy <- r
rd3dummy[] <- 0
rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
mice <- extract(ice, rd3points, method = "bilinear")
mice <- setValues(r, mice)
mice[is.na(mice)] <- 0
mice[mice > 99] <- NA
rm(rd3dummy, rd3points)

## Combine the masks
msk <- mask(mice, landmask)
msk[!is.na(msk)] <- 0

## Clean up
rm(landmask, ice, mice)

#------------------------------

# Function for intra-seasonal variance
intra_var <- function(raster_in) {
  raster_out <- raster::subset(raster_in, c(1:length(unique(dates$austral_year))))
  for (i in 1:length(unique(dates$austral_year))) {
    this_year <- dates$austral_year[i]
    dx <- which(dates$austral_year == this_year)
    this_stack <- raster::subset(raster_in, dx)
    raster_out[[i]]  <- raster::calc(this_stack, fun = sd, na.rm = T)
  }
  return(raster_out)
}

#-----------
# TOPOGRAPHY
#-----------
# Depth & Slope
DEPTH <- readbathy(topo = "gebco_19", xylim = extent(r_pad))
SLOPE <- terrain(DEPTH, opt = "slope")

DEPTH <- resample(DEPTH, r)
SLOPE <- resample(SLOPE, r)

DEPTH <- mask(DEPTH, msk)
SLOPE <- mask(SLOPE, msk)

# Distance to shelf
SHELFDIST <- readderivaadc("distance_shelf")
SHELFDIST <- resample(SHELFDIST, r)
SHELFDIST <- mask(SHELFDIST, msk)

# Distance to upper slope (Antarctica)
SLOPEDIST <- readderivaadc("distance_upper_slope")
SLOPEDIST <- resample(SLOPEDIST, r)
SLOPEDIST <- mask(SLOPEDIST, msk)

#-----------
# SST
#-----------
# SST
# High resolution SST (this fails, uses too much temp space)
if (FALSE) {
  # Note that GHRSST is only from 2002-06-01
  SST <- readghrsst(date = dates$dates[dates$dates > as.POSIXct("2002-06-01")], xylim = extent(r_pad))
  SST <- SST - 273.15 # K to C
}

# Or lower resolution SST (which allows fronts to be calculated)
if (TRUE) {
  SST <- readsst(date = dates$dates, xylim = extent(r_pad))
}

# SST front
SSTFRONT <- stack(SST)
for (i in 1:nlayers(SST)) {
  this_sst <- SST[[i]]
  SSTFRONT[[i]]  <- detectFronts(this_sst, method = "BelkinOReilly2009")
}

SSTFRONT <- mean(SSTFRONT)
SSTFRONT <- resample(SSTFRONT, r)
SSTFRONT <- mask(SSTFRONT, msk)

# SST intra-seasonal variance
SSTVAR <- intra_var(SST)
SSTVAR <- mean(SSTVAR)
SSTVAR <- resample(SSTVAR, r)
SSTVAR <- mask(SSTVAR, msk)

# Done with SST stack, calculate mean
SST <- mean(SST)
SST <- resample(SST, r)
SST <- mask(SST, msk)

#-----------
# EKE
#-----------
# Eddy kinetic energy
CURU <- readcurr(dates$dates, xylim = r, latest = FALSE, uonly = TRUE)
CURV <- readcurr(dates$dates, xylim = r, latest = FALSE, vonly = TRUE)

EKE <- 0.5*(CURU^2 + CURV^2)

EKE <- mean(EKE)
EKE <- resample(EKE, r)
EKE <- mask(EKE, msk)

#---------
# SSH
#---------
# SSH gradient
SSH <- readssh(date = dates$dates, xylim = extent(r_360_pad), ssha = F, lon180 = F)
SSH <- rotate(SSH)

SSHGRAD <- stack(SSH)
for (i in 1:nlayers(SSH)) {
  this_ssh <- SSH[[i]]
  SSHGRAD[[i]]  <- detectFronts(this_ssh, method = "BelkinOReilly2009")
}

SSHGRAD <- mean(SSHGRAD)
SSHGRAD <- resample(SSHGRAD, r)
SSHGRAD <- mask(SSHGRAD, msk)

# SSHA
SSHA <- readssh(date = dates$dates, xylim = extent(r_360_pad), ssha = T, lon180 = F)
SSHA <- rotate(SSHA)
SSHA <- mean(SSHA)
SSHA <- resample(SSHA, r)
SSHA <- mask(SSHA, msk)

#---------
# ICE
#---------
# Create a reference grid for calculating distance to ice edge
msk_ice <- readice_daily()
msk_ice[] <- 1
r_polar <- projectRaster(r, res = c(25000, 25000), crs = crs(msk_ice))
msk_wgs <- projectRaster(msk_ice, r, method = "ngb")

#-----------
# Concentration
ICE <- readice_daily(date = dates$dates)
ICE <- projectRaster(ICE, r)
ICE <- mask(ICE, msk_wgs, updatevalue = 0, updateNA = TRUE) # Set ice off the NSIDC grid to 0

# Distance
ICEDIST <- ICE
for (i in 1:length(dates$dates)) {
  print(paste(i, "of", length(dates$dates)), sep = " ")
  this_date <- dates$dates[i]
  this_ice <- readice_daily(date = this_date)
  edge <- rasterToContour(this_ice, levels = 15)
  
  # Distance to ice edge
  this_icedist <- r_polar
  # this_icedist <- this_ice
  # this_icedist[] <- NA
  dd <- rgeos::gDistance(edge, as(this_icedist, "SpatialPoints"), byid = TRUE)
  this_icedist[] = apply(dd,1,min)
  # Reproject both
  this_icedist <- projectRaster(this_icedist, r)
  this_ice <- projectRaster(this_ice, r)
  this_ice <- mask(this_ice, msk_wgs, updatevalue = 0, updateNA = TRUE)
  # Set distance within ice to negative
  values(this_ice)[which(values(this_ice) < 15)] <- 1
  values(this_ice)[which(values(this_ice) > 15)] <- -1
  this_icedist <- this_icedist * this_ice
  ICEDIST[[i]] <- this_icedist
}

ICEDIST <- mean(ICEDIST)
ICEDIST <- resample(ICEDIST, r)
ICEDIST <- mask(ICEDIST, msk)
ICEDIST <- ICEDIST/1000

# Intra-season variation
ICEVAR <- intra_var(ICE)

ICEVAR <- mean(ICEVAR)
ICEVAR <- resample(ICEVAR, r)
ICEVAR <- mask(ICEVAR, msk)
ICEVAR <- mask(ICEVAR, msk_wgs, updatevalue = 0, updateNA = TRUE)
ICEVAR <- mask(ICEVAR, msk)

ICE <- mean(ICE)
ICE <- resample(ICE, r)
ICE <- mask(ICE, msk)
ICE <- mask(ICE, msk_wgs, updatevalue = 0, updateNA = TRUE)
ICE <- mask(ICE, msk)

# Stack 'em
envars <- stack(DEPTH,
                SLOPE,
                SLOPEDIST,
                SHELFDIST,
                SST,
                SSTFRONT,
                SSTVAR,
                ICE,
                ICEDIST,
                ICEVAR,
                EKE,
                SSHA,
                SSHGRAD)

names(envars) <- c("DEPTH",
                   "SLOPE",
                   "SLOPEDIST",
                   "SHELFDIST",
                   "SST",
                   "SSTFRONT",
                   "SSTVAR",
                   "ICE",
                   "ICEDIST",
                   "ICEVAR",
                   "EKE",
                   "SSHA",
                   "SSHGRAD")

## Fill true data gaps using a focal filter
if (TRUE) {
  for (i in 1:nlayers(envars)) {
    foo <- envars[[i]]
    foo <- focal(foo, w=matrix(1,nrow=13,ncol=13), fun = "mean", na.rm = T)
    foo <- cover(SSTFRONT, foo)
    foo <- mask(foo, msk)
    envars[[i]] <- foo
  }
}

saveRDS(envars, "./output/envars.RDS")

# Create a dataframe for prediction
envar_frame_coords <- coordinates(envars)
envar_frame_envars <- as.data.frame(envars)
envar_frame <- cbind(envar_frame_coords, envar_frame_envars)
saveRDS(envar_frame, "./output/envar_frame.RDS")


#------------------------------
# Extract envars along real + simulated tracks

tracks <- readRDS("./output/tracks_real&amp;sim_region.RDS")
foo <- raster::extract(x = envars, y = as.data.frame(tracks[,c("lon", "lat")]))
tracks <- cbind(tracks, foo)
saveRDS(tracks, "./output/tracks_with_envars.RDS")