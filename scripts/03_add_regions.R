## Assign regions to real + simulated tracks

library(dplyr)

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

# Get tracks
tracks <- readRDS("./output/tracks_real&sim.RDS")

# Get metadata with region assigned, from regional variation work
meta <- read.csv("./data_in/metaWithBreedingStock&Region.csv",
                 stringsAsFactors = F)

# Filter
meta <- dplyr::filter(meta, meta$individual_id %in% unique(tracks$id))
meta$id <- meta$individual_id

# Combine
tracks <- left_join(x = tracks, y = dplyr::select(meta, id, region))

# Save output
saveRDS(tracks, "./output/tracks_real&sim_region.RDS")
