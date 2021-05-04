## Prepare and extract IWC whaling and survey data

library(dplyr)
library(data.table)
library(easycsv)
library(ggplot2)
library(SOmap)

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

#--------------------------------------------------------
## Whaling

# List the files
fls <- list.files(".\\data_in\\IWC\\IWCDBv6.1\\IndivData-CSVfmt\\", full.names = T)
fls <- fls[grepl(".csv", fls)]

# Read in all files separately
tst <- fread_folder(".\\data_in\\IWC\\IWCDBv6.1\\IndivData-CSVfmt\\")

# Remove North Atlantic and North Pacific
rm(`NA`, NP)

# Bind up the other datasets, note fill = T due to comments being
# read in as column names in some caes
dat_in <- rbind(IO, SA, SP, SHP1, SHP2, SU, fill = T)
rm(IO, SA, SP, SHP1, SHP2, SU, SHL)

dat_in <- dat_in[,1:38]

# Create a new dataframe, dealing with the empty and duplicated column names
dat <- data.frame(
  "species" = dat_in$Sp,
  "year" = dat_in$Yr,
  "month" = dat_in$Mon,
  "day" = dat_in$Day)
  
dat$lon_degrees <- dat_in$Lon
dat$lon_minutes <- dat_in[,21]$Mn
dat$lon_hemisphere <- dat_in$V22
dat$lon_accuracy <- dat_in[,23]$Ac

dat$lat_degrees <- dat_in$Lat
dat$lat_minutes <- dat_in[,17]$Mn
dat$lat_hemisphere <- dat_in$V18
dat$lat_accuracy <- dat_in[,19]$Ac

# Create coordinates
library(biogeo)
dat$lon <- dms2dd(dd = dat$lon_degrees, mm = dat$lon_minutes, ss = 0, ns = dat$lon_hemisphere)
dat$lat <- dms2dd(dd = dat$lat_degrees, mm = dat$lat_minutes, ss = 0, ns = dat$lat_hemisphere)

#-------------------------------
# Species
# These are coded as follows:
#   01 Pilot 	06 Sperm 	11 Right 	16 Pygmy Right 
# 02 Bottlenose 	07 Humpback 	12 Gray 	17 Cuvier's Beaked 
# 	03 Killer 	08 Sei 	13 Baird's Beaked 	18 Bowhead 
# 04 Blue 	09 Common Minke 	14 Baleen 	19 Beaked (unspecified)
# 05 Fin 	10 Bryde's 	15 Pygmy Blue 	20 Antarctic Minke

#-------------------------------
# Accuracy:
# The fields denoted 'Ac' following the catch position contain an indicator defining how accurately the position was reported.
# 0: Unknown position
# 1: Exact position given to nearest minute
# 2: Exact position given to nearest degree (or half degree)
# 3: Approximate position
# 4: position calculated from distance and bearing

#-------------------------------
# Filter

# Species
dat <- filter(dat, species == 7)

# Accuracy
if (TRUE) {
  # Exact positions to nearest degree
  dat <- filter(dat, lon_accuracy == 2 | lon_accuracy == 1)
  dat <- filter(dat, lat_accuracy == 2 | lat_accuracy == 1)
} else {
  # Or only the most accurate
  dat <- filter(dat, lat_accuracy == 1 & lon_accuracy == 1)
}

# South of 40S
dat <- filter(dat, lat <= -40)

# Only in given months
which_months <- c(11, 12, 1, 2, 3)
dat <- filter(dat, month %in% which_months)

#-------------------------------
# Plot
ggplot(data = dat, aes(x = lon, y = lat)) +
  # geom_point(alpha = 0.1) +
  coord_quickmap() +
  geom_hex() +
  scale_fill_viridis_c()

# In SOmap
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(dat$lon, dat$lat, pch = 16)

#-------------------------------
# Binned on a grid
library(raster)
library(pals)

# Create a reference grid
grd <- raster(ext = extent(-180, +180, -80, -40), res = c(1, 1), crs = "+proj=longlat +datum=WGS84 +no_defs")

# Create an sp class object from the whaling data
dat_sp <- dat
coordinates(dat_sp) <- dat[ , c("lon", "lat")]
crs(dat_sp) <- "+proj=longlat +datum=WGS84 +no_defs"

# Count on the grid
dat_grid <- rasterize(dat_sp, grd, fun='count')[["ID"]]

# Plot
tiff("./figures/mapWhalingCatches.tiff",
     height = 5,
     width = 5,
     units = "in",
     res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

# SOplot(dat_grid, col = ocean.thermal(125))
SOplot(dat$lon, dat$lat, pch = 16, cex = 0.3, col = "#EE7733")

dev.off()

# Write raster
dat_grid[!is.na(dat_grid) & dat_grid > 1] <- 1
writeRaster(dat_grid, "./output/iwc_catches.grd", overwrite = T)

#--------------------------------------------------------
## Surveys

surv <- read.csv("./data_in/IWC/Allsght.csv", stringsAsFactors = F)

# Filter
# Humpacks = species code 7
surv <- dplyr::filter(surv, SpeciesCode1 == 7)
#Months
surv <- filter(surv, Month %in% which_months)

# Fix year
surv[which(surv$Year > 50), ]$Year <- surv[which(surv$Year > 50), ]$Year + 1900
surv[which(surv$Year < 50), ]$Year <- surv[which(surv$Year < 50), ]$Year + 2000


# Create coordinates
# Assume columns 'DecimalLongitude' and 'DecimalLatitude' are the decimal portion of minutes
surv[is.na(surv$DecimalLatitudr), "DecimalLatitudr"] <- 0 # if missing set to 0
surv[is.na(surv$DecimalLongitude), "DecimalLongitude"] <- 0 # if missing set to 0

surv$MinutesLatitude <- as.numeric(paste(surv$MinutesLatitude, surv$DecimalLatitudr, sep = ".")) # Add decimal minutes comp.
surv$MinutesLongitude <- as.numeric(paste(surv$MinutesLongitude, surv$DecimalLongitude, sep = ".")) # Add decimal minutes comp.

surv$lat <- dms2dd(dd = surv$DegreesLatitude, mm = surv$MinutesLatitude, ss = 0, ns = surv$LatitudeAttribute)
surv$lon <- dms2dd(dd = surv$DegreesLongitude, mm = surv$MinutesLongitude, ss = 0, ns = surv$LongitudeAttribute)

# Filter to -40
surv <- surv[surv$lat <= -40, ]

# Plot to check
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)
SOplot(surv$lon, surv$lat, pch = 16)

# Create an sp class object from the survey data
surv_sp <- surv
coordinates(surv_sp) <- surv[ , c("lon", "lat")]
crs(surv_sp) <- "+proj=longlat +datum=WGS84 +no_defs"

# Count on the grid
surv_grid <- rasterize(surv_sp, grd, fun='count')[["ID"]]

# Plot
tiff("./figures/mapWhalingSurvey.tiff",
     height = 5,
     width = 5,
     units = "in",
     res = 300)
SOmap(trim = -40,
      bathy_legend = FALSE,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

# SOplot(surv_grid, col = ocean.thermal(125))
SOplot(surv$lon, surv$lat, pch = 16, cex = 0.3, col = "#EE7733")

dev.off()

# Write raster
surv_grid[!is.na(surv_grid) & surv_grid > 1] <- 1
writeRaster(surv_grid, "./output/iwc_survey.grd", overwrite = T)


# Combine catches and survey
dat_grid[is.na(dat_grid)] <- 0
surv_grid[is.na(surv_grid)] <- 0
iwc_grid <- overlay(surv_grid, dat_grid, fun = sum)

# Write raster
iwc_grid[!is.na(iwc_grid) & iwc_grid > 1] <- 1
writeRaster(iwc_grid, "./output/iwc_catches_and_survey.grd", overwrite = T)

#------------------------------------------------------------
# Alternative maps
prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Catch data
ll_catch <- dat
coordinates(ll_catch) <- c("lon", "lat")
projection(ll_catch) <- "+proj=longlat +datum=WGS84"
points_catch <- SOproj(ll_catch, target = prj)

# Survey data
ll_survey <- surv
coordinates(ll_survey) <- c("lon", "lat")
projection(ll_survey) <- "+proj=longlat +datum=WGS84"
points_survey <- SOproj(ll_survey, target = prj)

# Longtidue labels
lon_labels <- data.frame("lon" = c(seq(-180, +135, 45)),
                         "lat" = -34,
                         "lon_name" = c("180°W|180°E", "135°W", "90°W", "45°W", "0°W|0°E", "45°E", "90°E", "135°E"))
coordinates(lon_labels) <- c("lon", "lat")
projection(lon_labels) <- "+proj=longlat +datum=WGS84"
lon_labels <- SOproj(lon_labels, target = prj)


this_map <- SOgg(SOmap(trim = -40,
                       bathy_legend = FALSE,
                       border_col = c("white", "white"),
                       border_width = 0.01,
                       straight = TRUE,
                       graticules = TRUE))


p_a <- plot(this_map) +
  geom_point(data = as.data.frame(points_catch), aes(lon, lat),
             size = 0.7, colour = "#EE7733") +
  geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name),
            colour = "black", size = 3) +
  guides(fill = "none") +
  labs(subtitle = "a")

p_b <- plot(this_map) +
  geom_point(data = as.data.frame(points_survey), aes(lon, lat),
             size = 0.7, colour = "#EE7733") +
  geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name),
            colour = "black", size = 3) +
  guides(fill = "none") +
  labs(subtitle = "b")


tiff("./figures/iwc_catch_map.tiff",
     height = 5,
     width = 5.5,
     units = "in",
     res = 300)
print(p_a)
dev.off()

tiff("./figures/iwc_survey_map.tiff",
     height = 5,
     width = 5.5,
     units = "in",
     res = 300)
print(p_b)
dev.off()
