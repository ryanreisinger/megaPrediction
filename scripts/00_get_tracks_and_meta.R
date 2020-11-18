# Copy combined tracking data and metadata into working folder
# Ryan Reisinger
# September 2020

library(dplyr)

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

# -------------
## Tracks
# -------------

trk.fls <- list.files("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/mega/data_formatted/tracks/", full.names = T)
trk.fls <- trk.fls[grepl(".csv", trk.fls)]

tracks <- do.call(rbind, lapply(trk.fls, read.csv, stringsAsFactors = F))

# -------------
# Replace old location classes
tracks <- tracks %>%
  mutate(location_quality = ifelse(location_quality == "-9", "Z", location_quality)) %>%
  mutate(location_quality = ifelse(location_quality == "-3", "Z", location_quality)) %>%
  mutate(location_quality = ifelse(location_quality == "-2", "B", location_quality)) %>%
  mutate(location_quality = ifelse(location_quality == "-1", "A", location_quality)) %>%
  mutate(location_quality = ifelse(location_quality == "Z", "B", location_quality))

# -------------
# Filter individuals not to use
tracks <- dplyr::filter(tracks, individual_id != 112694)
tracks <- dplyr::filter(tracks, individual_id != 20683)
tracks <- dplyr::filter(tracks, individual_id != "Mn_WAVES14_Lander01")
tracks <- dplyr::filter(tracks, individual_id != "Entangled whale")

# -------------
# Drop the handful of records with no date information
tracks <- filter(tracks, !is.na(tracks$date))

# -------------
# Drop records in the far northern hemisphere
# These are all test locations or errors
tracks <- dplyr::filter(tracks, tracks$decimal_latitude < 20)

# -------------
# In Brazil data, drop 'Tagging' records
tracks <- dplyr::filter(tracks, tracks$location_quality != "Tagging")

# -------------
# Write output
write.csv(tracks, "./data_in/tracks.csv", row.names = F)

# -------------
## Metadata
# -------------

met.fls <- list.files("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/mega/data_formatted/meta/", full.names = T)
met.fls <- met.fls[grepl(".csv", met.fls)]

meta <- do.call(rbind, lapply(met.fls, read.csv, stringsAsFactors = F))

# Add entries for tracks that have no metadata
meta$meta_source <- "meta" # add a flag indicating if metadata were generated from tracks

`%nin%` = Negate(`%in%`) # Negation

foo.tracks <- unique(tracks$individual_id) # Track ids
foo.meta <- unique(meta$individual_id) # Meta ids
foo.new <- foo.tracks[which(foo.tracks %nin% foo.meta)] # Which track ids are not in meta
foo.new <- foo.new[! is.na(foo.new)] # Drop NAs

# Make extra metadata id by id
extra.meta <- data.frame()

for (i in 1:length(foo.new)) {
  this.id <- foo.new[i]
  this.dat <- filter(tracks, tracks$individual_id == this.id)
  this.dat <- this.dat[1, ]
  this.meta <- data.frame("dataset_identifier" = this.dat$dataset_identifier,
                          "data_owner" = NA,
                          "contact_email" = NA,
                          "file_name" = NA,
                          "individual_id" = this.dat$individual_id,
                          "device_id" = this.dat$device_id,
                          "device_type" = this.dat$device_type,
                          "year" = NA,
                          "month" = NA,
                          "day" = NA,
                          "time" = NA,
                          "time_zone" = NA,
                          "deployment_site" = NA,
                          "deployment_decimal_latitude" = NA,
                          "deployment_decimal_longitude" = NA,
                          "sex" = NA,
                          "how_sexed" = NA,
                          "age_class" = NA,
                          "genotyped" = NA,
                          "age." = NA,
                          "progesterone" = NA,
                          "If.yes..status" = NA,
                          "comments" = NA,
                          "meta_source" = "track")
  extra.meta <- bind_rows(extra.meta, this.meta)
}

meta <- bind_rows(meta, extra.meta)

# Drop metadata for which there is no tracking data
meta <- filter(meta, meta$individual_id %in% foo.tracks)

# Write a quick combined copy
write.csv(meta, "./data_in/meta.csv", row.names = F)
