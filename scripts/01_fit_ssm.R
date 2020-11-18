# Fit State-Space Model to Tracking Data
# Ryan Reisinger

library(dplyr)
library(tidyr)
library(foieGras)
library(sf)
library(SOmap)

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

## Get data
tracks <- read.csv("./data_in/tracks.csv", stringsAsFactors = FALSE)

# -------------
## Organise for foieGras

## Create a dataframe for foieGras
dat <- dplyr::select(tracks,
                     individual_id,
                     date,
                     location_quality,
                     decimal_longitude,
                     decimal_latitude,)

dat <- dplyr::rename(dat,
                     id = individual_id,
                     date = date,
                     lc = location_quality,
                     lon = decimal_longitude,
                     lat = decimal_latitude)

## Remove ids with less than 10 locations
dat <- dplyr::group_by(dat, id) %>% 
  filter(n() > 9L) %>% 
  ungroup(.)

## Prefilter individually to avoid errors with projection
ids <- unique(dat$id)
all.d <- data.frame()

for (i in ids) {
  print(i)
  this.d <- dat[dat$id == i, ]
  this.d <- this.d[complete.cases(this.d), ]
  this.d <- fit_ssm(this.d, vmax = 5, pf = TRUE)
  this.d <- st_transform(this.d, crs = "+proj=longlat +datum=WGS84")
  pts <- st_coordinates(this.d)
  that.d <- as.data.frame(this.d)[ , c("id", "date", "lc", "keep")]
  that.d$lon <- pts[,1]
  that.d$lat <- pts[,2]
  that.d <- that.d[that.d$keep == "TRUE", ]
  that.d$keep <- NULL
  all.d <- rbind(all.d, that.d)
}

dat <- all.d

#--------------------------------------------------------------------------
## Fit the SSM in foieGras

  
  ## Fit the model
  fit <- foieGras::fit_ssm(dat, pf = FALSE, spdf = FALSE, model = "rw", time.step = 24)
  
  ## Grab the output
  out <- grab(fit, "predicted", as_sf = FALSE)
  
  ## Filter below 40 degrees
  out <- dplyr::filter(out, lat < -40)
  
  # Check a plot of each track
  ids <- unique(out$id)
  
  ## Plot
  pdf("./output/checkTracksFoieGras.pdf", paper = "a4", useDingbats = F)
  for (i in ids) {
    print(i)
    this.d <- out[out$id == i, ]
    print(SOmap_auto(this.d$lon, this.d$lat), main = i)
  }
  dev.off()

saveRDS(out, "./output/foieGras_ssm.RDS")