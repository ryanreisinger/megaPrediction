## Simulate pseudotracks

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

# install.packages("remotes") # Install 'remotes' if you don't have it
# library(remotes)
# install_github("AustralianAntarcticDataCentre/availability") # install 'availability'
library(availability)
library(dplyr)

## Get fitted tracks
tracks <- readRDS("./output/foieGras_ssm.RDS")

#------------------------------------------
# Some helper functions

# Simple helper to plot tracks
plotTracks <- function(..., pal = c("black", "dodgerblue")) {
  require(maptools)
  data(wrld_simpl)
  trks <- list(...)
  xlim <- do.call(range, lapply(trks, function(trk) trk[, 1]))
  ylim <- do.call(range, lapply(trks, function(trk) trk[, 2]))
  plot(xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xlab = "longitude", ylab = "latitude")
  for (k in seq(from = floor((min(xlim)+180)/360), to = floor((max(xlim)+180)/360), by = 1))
    plot(elide(wrld_simpl, shift = c(k*360, 0)), col = "grey30", border = "grey30",
         add = TRUE, xlim = xlim, ylim = ylim)
  for(k in seq_along(trks)) points(trks[[k]], col = pal[k])
}

# A custom landmask function
# WARNING: this approach is VERY slow!
# Eventual aim is to check sea ice here
require(raster)
r <- raster("./data_landmask/ETOPO1_custom.tif")

myMask <- function() {
  function(tm, pt) {
    xpt <- pt[1]
    ypt <- pt[2]
  v <- raster::extract(r, matrix(c(xpt, ypt), nrow = 1, ncol = 2))
  if (v > 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
  }
}

#------------------------------------------
# Create a function that does multiple simulations
multi_sims <- function(track, nsims = 50) {
  hold <- data.frame()
  for (j in 1:nsims) {
    cat("Simulating", j, "of", nsims, "\n")
    this_arfit <- surrogateARModel(dplyr::select(track, lon, lat))
    this_sim <- surrogateAR(model = this_arfit, # This is the fitted model
                            xs = dplyr::select(track, lon, lat), # This is the 'template' track
                            fixed = track$fixed, # These are the fixed points
                            # The built-in land mask:
                            # point.check = gshhsMask(res = 0.05), # This is the land mask
                            # Or a custom landmask relying on a custom geotiff:
                            point.check = etopoMask(basename = "ETOPO1_custom", path = "./data_landmask", tmp = "./data_landmask",
                                                    land = FALSE),
                            # Or a custom written function (very slow!):
                            # point.check = myMask(),
                            partial = TRUE # Should a partial track be allowed
    )
    this_sim <- data.frame("id" = unique(track$id),
                           "lon" = this_sim$xs[,1],
                           "lat" = this_sim$xs[,2],
                           "date" = track$date,
                           "n_sim" = j)
    if (nrow(this_sim) > 0) {
      hold <- rbind(hold, this_sim)
    }
  }
  return(hold)
}

#------------------------------------------
# Loop over each individual, applying this function to each

# Try filtering by track length
tracks <- dplyr::group_by(tracks, id) %>% 
  filter(n() > 5L) %>% 
  ungroup(.)

# Set the parameters
nsims = 50

# Define the IDs
ids <- unique(tracks$id)

# A dataframe to hold the output
master <- data.frame()

# Then loop through each id
for (i in 1:length(ids)) {
  
  # Individual i
  this_id <- ids[i]
  print(this_id)
  
  # Get track for individual i
  this_track <- dplyr::filter(tracks, id == this_id)
  
  # Set fixed location
  this_track$fixed <- FALSE
  this_track$fixed[1] <- TRUE
  
  # Simulate
  sims <- multi_sims(track = this_track, nsims = nsims)
  
  # Add to previous sims if all did not fail
  if (nrow(sims) > 0) {
    master <- rbind(master, sims)
  } else {
    cat("All sims failed")
  }
  
}

# Plot
plotTracks(dplyr::select(master, lon, lat), # The original tracks
           dplyr::select(tracks, lon, lat), # The simulated tracks
           pal = c("grey", "dodgerblue") # Colour palette
)

# Combine real + simulated tracks
tracks$n_sim <- 0
tracks$real <- "O"
master$real <- "S"

track_out <- rbind(dplyr::select(tracks, id, date, lon, lat, n_sim, real),
                   dplyr::select(master, id, date, lon, lat, n_sim, real))

track_out$track_id <- paste(track_out$id, track_out$n_sim, sep = "__")

# Save output
saveRDS(track_out, "./output/tracks_real&sim.RDS")
