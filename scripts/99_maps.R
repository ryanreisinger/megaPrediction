# Maps

library(raster)
library(SOmap)
library(pals)
library(ggplot2)
library(ggnewscale)

setwd("C:\\Users\\Ryan Reisinger\\Documents\\Academic\\UCSC\\Work\\Analysis\\megaPrediction")

#-------------------------------------
# Setup
#-------------------------------------

#-------------------------------------
# Tracks & covariates
#-------------------------------------

# Get data
envars <- readRDS("./output/envars.RDS")
tracks <- readRDS("./output/tracks_with_envars.RDS")

#-------------------------------------
# Plot stuff
#-------------------------------------
prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

lon_labels <- data.frame("lon" = c(seq(-180, +135, 45)),
                         "lat" = -33,
                         "lon_name" = c("180°W|180°E", "135°W", "90°W", "45°W", "0°W|0°E", "45°E", "90°E", "135°E"))
coordinates(lon_labels) <- c("lon", "lat")
projection(lon_labels) <- "+proj=longlat +datum=WGS84"
lon_labels <- SOproj(lon_labels, target = prj)

#-------------------------------------
# Environmental covariates
#-------------------------------------

envar_titles <- data.frame("name" = names(envars),
                           "title" = c("Depth\n(m)",
                                       "Bottom slope\n(deg)",
                                       "Distance to slope\n(km)",
                                       "Distance to shelf\n(km)",
                                       "Sea surface temperature\n(°C)",
                                       "Sea surface temperature gradient\n(°C/km)",
                                       "Sea surface temperature variance\n(°C)",
                                       "Sea ice concentration\n(%)",
                                       "Distance to sea ice edge\n(km)",
                                       "Sea ice concentration variance\n(%)",
                                       "Eddy kinetic energy\n(unit)",
                                       "Sea surface height anomaly\n(m)",
                                       "Sea surface height gradient\n(unit)"))

# Option 1
# SOplot
this_map <- SOmap(trim = -40,
                  bathy_legend = FALSE,
                  border_col = c("white", "white"),
                  border_width = 0.01,
                  straight = TRUE,
                  graticules = TRUE)

# Get rid of bathy
raster::values(this_map$bathy[[1]]$plotargs$x) <- NA_real_
this_map$bathy_legend <- NULL

# Plot
for (i in 1:nlayers(envars)) {
tiff(paste0("./figures/envars/envars_", names(envars[[i]]), ".tiff"),
     height = 2.5,
     width = 5,
     units = "in",
     res = 300)
plot(this_map)
SOplot(SOproj(envars[[i]]),
       col = parula(125),
       legend.args = list(text = envar_titles$title[i]))
dev.off()
}

# Option 2
for (i in 1:nlayers(envars)) {
  # i = 1
  this_raster <- envars[[i]]
  this_var <- names(envars[[i]])
  this_caption <- envar_titles$title[i]
  names(this_raster) <- c("val")
  
  this_map <- SOgg(SOmap(trim = -40,
                         bathy_legend = FALSE,
                         border_col = c("white", "white"),
                         border_width = 0.01,
                         straight = TRUE,
                         graticules = TRUE))
  
  this_map$bathy <- NULL

  
  p <- plot(this_map) +
    new_scale_fill() +
    geom_raster(data = as.data.frame(SOproj(this_raster), xy = TRUE),
                aes(x = x, y = y, fill = val)) +
    scale_fill_gradientn(colors = ocean.haline(125), na.value = NA, name = this_caption) +
    labs(subtitle = this_var) +
    geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name), colour = "black", size = 2) +
    theme(legend.position = "bottom") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 10))
  
  tiff(paste0("./figures/envars_alt/envars_alt_", this_var, ".tiff"),
       height = 5,
       width = 5,
       units = "in",
       res = 600)
  print(p)
  dev.off()
}

#-------------------------------------
# Predictions
#-------------------------------------

these_models <- c("Mr_Atlantic", "Mr_EastIndian", "Mr_EastPacific", "Mr_Pacific", "Mr_WestPacific",
                  "M1", "M2", "M3", "M4", "M5")

for (i in 1:length(these_models)) {

# which_model <- "M1"
  which_model <- these_models[i]

#--------------------------

# Get the raster
this_file <- paste0("./output/predictions/", which_model, "_raster.grd")
r <- raster(this_file)
crs(r) <- "+init=epsg:4326"

names(r) <- c("val")

p <- plot(this_map) +
  new_scale_fill() +
  geom_raster(data = as.data.frame(SOproj(r), xy = TRUE),
              aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(colors = parula(125), na.value = NA, name = "p(Observed track)",
                       limits = c(0, 1)) +
  geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name), colour = "black", size = 4) +
  labs(subtitle = which_model) +
  theme(plot.subtitle = element_text(size = 13))

tiff(paste0("./figures/predictions/predictions_", which_model, ".tiff"),
     height = 5,
     width = 6.5,
     units = "in",
     res = 600)
print(p)
dev.off()

}

#-------------------------------------
# Tracks
#-------------------------------------
tracks <- readRDS("./output/tracks_with_envars.RDS")

# tracks <- tracks[tracks$real == "O", ]
# 
# ll_tracks <- tracks
# coordinates(ll_tracks) <- c("lon", "lat")
# projection(ll_tracks) <- "+proj=longlat +datum=WGS84"
# points_tracks <- SOproj(ll_tracks, target = prj)

this_map <- SOgg(SOmap(trim = -40,
                       bathy_legend = FALSE,
                       border_col = c("white", "white"),
                       border_width = 0.01,
                       straight = TRUE,
                       graticules = TRUE))

# Real + simulated
ll_tracks <- tracks
coordinates(ll_tracks) <- c("lon", "lat")
projection(ll_tracks) <- "+proj=longlat +datum=WGS84"
points_tracks <- SOproj(ll_tracks, target = prj)

p_a <- plot(this_map) +
  geom_point(data = as.data.frame(points_tracks[points_tracks$real == "S", ]), aes(x = lon, y = lat),
             size = 0.7, alpha = 1, col = "grey") +
  geom_point(data = as.data.frame(points_tracks[points_tracks$real == "O", ]), aes(x = lon, y = lat),
             size = 0.7, alpha = 1, col = "black") +
  geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name),
            colour = "black", size = 2) +
  # scale_colour_manual(values = c("black", "grey"),
  #                     name = "Track",
  #                     labels = c("Observed", "Simulated")) +
  guides(colour = guide_legend(override.aes = list(size = 3)), fill = FALSE) +
  
  guides(fill = "none") +
  labs(subtitle = "a") +
  theme(plot.subtitle = element_text(size = 13))

tiff("./figures/tracks_all.tiff",
     height = 5,
     width = 6.5,
     units = "in",
     res = 600)
print(p_a)
dev.off()


# Real only, by region
ll_tracks <- tracks

#Filter
ll_tracks <- ll_tracks[ll_tracks$real == "O", ]
ll_tracks <- ll_tracks[!is.na(ll_tracks$region), ]

coordinates(ll_tracks) <- c("lon", "lat")
projection(ll_tracks) <- "+proj=longlat +datum=WGS84"
points_tracks <- SOproj(ll_tracks, target = prj)

p_b <- plot(this_map) +
  geom_point(data = as.data.frame(points_tracks), aes(x = lon, y = lat, col = region),
             size = 0.7) +
  geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name),
            colour = "black", size = 2) +
  scale_colour_manual(values = brewer.set1(5),
                      name = "Region",
                      labels = c("Atlantic", "East Indian", "East Pacific", "Central Pacific", "West Pacific")) +
  guides(colour = guide_legend(override.aes = list(size = 3)), fill = FALSE) +
  
  guides(fill = "none")

tiff("./figures/tracks_regions.tiff",
     height = 5,
     width = 6.5,
     units = "in",
     res = 600)
print(p_b)
dev.off()

# Real only, individual maps by region
regions <- unique(tracks$region)
regions <- regions[!is.na(regions)]

for (i in 1:length(regions)) {
  this_region <- regions[i]
  this_data <- as.data.frame(points_tracks)
  this_data <- this_data[this_data$region == this_region, ]
  
  if(this_region == "Atlantic") {
  this_colour <- "#e41a1c"
  }
  if(this_region == "EastIndian") {
    this_colour <- "#377eb8"
  }
  if(this_region == "WestPacific") {
    this_colour <- "#ff7f00"
  }
  if(this_region == "Pacific") {
    this_colour <- "#984ea3"
  }
  if(this_region == "EastPacific") {
    this_colour <- "#4daf4a"
  }
  
  
  p <- plot(this_map) +
    geom_point(data = this_data, aes(x = lon, y = lat),
               size = 0.7,
               col = this_colour) +
    geom_text(data = as.data.frame(lon_labels), aes(x = lon, y = lat, label = lon_name),
              colour = "black", size = 4) +
        guides(fill = "none") +
    labs(subtitle = this_region) +
    theme(plot.subtitle = element_text(size = 13))
  
  tiff(paste0("./figures/tracks_regions_", this_region, ".tiff"),
       height = 5,
       width = 6.5,
       units = "in",
       res = 600)
  print(p)
  dev.off()
  
}
