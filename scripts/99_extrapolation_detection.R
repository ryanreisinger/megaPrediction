# Extrapolation detection for regional models

library(dsmextra)

setwd("/mnt/home/ryan/humpbacks/megaPrediction")

#---------------------------------
# Reference data
# the calibration space
refdat <- readRDS("./output/tracks_with_envars.RDS")

refdat <- refdat[,c(3, 4, 8:18)]
refdat <- refdat[complete.cases(refdat), ]

#---------------------------------
# Projection data
# the space where predictions are made
projdat <- readRDS("./output/envar_frame.RDS")

projdat$EKE <- NULL
projdat$SSHA <- NULL
projdat$SSHGRAD <- NULL

# projdat <- projdat[complete.cases(projdat),]

#---------------------------------
# Names
envars <- names(projdat)[-c(1:2)]
proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#---------------------------------
# Function to check extrapolation for a given region
ex_foo <- function(region = "Atlantic") {
ex <- compute_extrapolation(samples = refdat[refdat$region == region, ],
                     prediction.grid = projdat,
                     covariate.names = envars,
                     coordinate.system = proj,
                     resolution = 0.1)

# Total cells
tot <- length(ex$data$all$ExDet)

# Combinatorial extrapolation
com_ex <- length(ex$data$all$ExDet[ex$data$all$ExDet > 1]) / tot
com_ex <- com_ex * 100

# Univariate extrapolation
uni_ex <- length(ex$data$all$ExDet[ex$data$all$ExDet < 0]) / tot
uni_ex <- uni_ex * 100

# # Univariate extrapolation %
# ex$summary$extrapolation$univariate.p

foo <- data.frame("region" = region,
                  "Univariate_ex_p" = uni_ex,
                  "Combinatorial_ex_p" = com_ex)

return(foo)

}

#---------------------------------
# Calculate each
these_regions <- unique(refdat$region)

results <- do.call(rbind, lapply(these_regions, ex_foo))

# Append to model table
mod_results <- read.csv("./output/model_results_table.csv", stringsAsFactors = F)
results$Model <- paste0("Mr_", results$region)
results$region <- NULL

results <- merge.data.frame(x = mod_results, y = results, all.x = T)

# Calculate correlations
foo <- results[!is.na(results$Combinatorial_ex_p), ]
cor(foo$Internal_cross_validation_mean, foo$Univariate_ex_p, method = "spearman")
cor(foo$Internal_cross_validation_mean, foo$Combinatorial_ex_p, method = "spearman")


# Tidy up
results$Combinatorial_ex_p <- round(results$Combinatorial_ex_p, 2)
results$Univariate_ex_p <- round(results$Univariate_ex_p, 2)

# Write output
write.csv(results, "./output/model_results_table_w_exdet.csv", row.names = F)