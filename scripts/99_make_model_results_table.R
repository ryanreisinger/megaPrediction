# Create model results table

setwd("C:\\Users\\Ryan Reisinger\\Documents\\Academic\\UCSC\\Work\\Analysis\\megaPrediction")

library(dplyr)

# Get results
fls <- list.files("./output/fitted_models_results/", full.names = T)
dat <- do.call(rbind, lapply(fls, read.csv, stringsAsFactors = F))

# Organise
out <- data.frame("Model" = dat$model,
                  "N_tracks" = dat$n_tracks,
                  "Internal_cross_validation_mean" = dat$ROC,
                  "Internal_cross_validation_sd" = dat$ROCSD,
                  "External_validation_all_tracks" = dat$validation_all_tracks,
                  "External_validation_catches_and_sightings" = dat$validation_iwc,
                  "Rank" = rank(-dat$validation_iwc))

# Tidy up
out$Internal_cross_validation_mean <- round(out$Internal_cross_validation_mean, 3)
out$Internal_cross_validation_sd <- round(out$Internal_cross_validation_sd, 3)
out$External_validation_all_tracks <- round(out$External_validation_all_tracks, 3)
out$External_validation_catches_and_sightings <- round(out$External_validation_catches_and_sightings, 3)

# Save
write.csv(out, "./output/model_results_table.csv", row.names = F)

# Correlations
dat_r <- dplyr::filter(out,
                       Model == "Mr_Atlantic" |
                         Model == "Mr_EastIndian" |
                         Model == "Mr_EastPacific" |
                         Model == "Mr_Pacific" |
                         Model == "Mr_WestPacific")

cor(dat_r$N_tracks, dat_r$External_validation_all_tracks, method = "spearman")
cor(dat_r$N_tracks, dat_r$External_validation_catches_and_sightings, method = "spearman")
