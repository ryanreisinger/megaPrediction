## Fit models

setwd("/mnt/home/ryan/humpbacks/megaPrediction")

## Various machine learning models fit in caret
library(ranger)
library(caret)
library(e1071)
library(pROC)

# Utils
library(plyr)
library(tidyr)
library(ggplot2)

## Libraries for parallel processing
library(iterators)
library(foreach)
library(parallel)
library(doParallel)

# Raster
library(raster)

#-------------------------------------------------------------------
# Load tracking data
dat <- readRDS("./output/tracks_with_envars.RDS")

# For caret, the response must be factor
dat$real <- as.factor(dat$real)

# Depending on the model, caret can only run with complete cases,
# this also means EKE, SSHA and SSHGRAD are excluded
dat <- dat[,1:19]
dat <- dat[complete.cases(dat[,9:19]), ]

#-------------------------------------------------------------------
# Load data for prediction
area_grid <- readRDS("./output/envar_frame.RDS")
area_grid$EKE <- NULL
area_grid$SSHA <- NULL
area_grid$SSHGRAD <- NULL

# Add IWC data for validation
iwc <- raster("./output/iwc_catches_and_survey.grd")
area_grid$iwc <- raster::extract(x = iwc, y = area_grid[,c("x", "y")])
area_grid$iwc[which(area_grid$iwc == 1)] <- "O"
area_grid$iwc[which(area_grid$iwc == 0)] <- "S"

#-------------------------------------------------------------------
# Modelling
#-------------------------------------------------------------------

## General parameters common across all models

allow_par <- FALSE
number_k <- 10
number_repeats <- 30
number_trees <- 2000

# Parameter grid
rfGrid <-  data.frame(mtry = c(2, 3, 4),
                      splitrule = "gini",
                      min.node.size = 1)

# Figure params
varimp_fig_height = 60/25.4
varimp_fig_width = 75/25.4

# Function for varimp
varimpR <- function(mod = M0, mod_name = "M0") {
tmp <- varImp(mod)$importance
tmp$Covariate <- row.names(tmp)
tmp$Model <- mod_name
tmp <- tmp[order(tmp$Overall, decreasing = TRUE), ]
write.csv(tmp, paste0("./output/fitted_models_varimp/varimp_", mod_name, ".csv"), row.names = F)

tmp$Covariate <- factor(tmp$Covariate, levels = rev(tmp$Covariate))
p <- ggplot(data = tmp) +
  geom_segment(aes(x = Covariate, xend = Covariate, y=0, yend = Overall), colour = "black") +
  geom_point(aes(x = Covariate, y = Overall), colour = "black") +
  coord_flip() +
  labs(x = "", y = "Relative variable importance", subtitle = unique(tmp$Model)) +
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"))

pdf(paste0("./figures/varimp/", mod_name, ".pdf"), height = varimp_fig_height, width = varimp_fig_width, useDingbats = FALSE)
print(p)
dev.off()

rm(tmp, p)
}


#--------------------------------------------
# 1. The naive global model - M1
#--------------------------------------------
# Here we use all data, fit a single model, and predict to the whole area

# Select our data
d1 <- dat # Using everything

# Create folds
folds <- groupKFold(group = d1$id, k = number_k) # This is caret's new built-in function

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)

# Set up parallel
# clust <- makeCluster(detectCores() - 1) # leave 1 core for OS
# registerDoParallel(clust)

# Fit the model
system.time(
  M1 <- train(x = as.data.frame(dplyr::select(d1,
                                              DEPTH,
                                              SLOPE,
                                              SLOPEDIST,
                                              SHELFDIST,
                                              SST,
                                              SSTFRONT,
                                              SSTVAR,
                                              ICE,
                                              ICEDIST,
                                              ICEVAR)),
              y = d1$real,
              method = "ranger",
              metric = "ROC",
              trControl = tc,
              tuneGrid = rfGrid,
              importance = "impurity",
              num.trees = number_trees)
)

# Inspect
M1

# Variable importance
varImp(M1)

# Run function to plot and save varimp
varimpR(mod = M1, mod_name = "M1")

# Save
saveRDS(M1, "./output/fitted_models/M1.RDS")

# Stop the cluster
# stopCluster(clust)
# registerDoSEQ()

#-----------------
# Predict
# Performance
dat$M1_predicted_prob <- predict(M1, newdata = dat, type = "prob")$O
dat$M1_predicted_class <- predict(M1, newdata = dat, type = "raw")

internal_performance <- M1$results[which(M1$results$ROC == max(M1$results$ROC)),]$ROC
internal_performance_sd <- M1$results[which(M1$results$ROC == max(M1$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_p <- NA
area_grid$M1_p[dx] <- predict(M1, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_p")])
writeRaster(r, "./output/predictions/M1_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M1$results[which(M1$results$ROC == max(M1$results$ROC)),]
results$model <- "M1"
results$n_tracks <- length(unique(d1$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1.csv", row.names = F)
rm(results)

#--------------------------------------------
# 2. Unweighted ensemble / simple mean - M2
#--------------------------------------------
# Here we fit models to each region, and for the global area take the mean of their predictions

# First, we fit the regional models
#----------------------
# a. Pacific

# Data
d2_Pacific <- dplyr::filter(dat, region == "Pacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
         DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR)

# Folds
folds <- groupKFold(group = d2_Pacific$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
Mr_Pacific <- train(x = as.data.frame(dplyr::select(d2_Pacific,
                                                    DEPTH,
                                                    SLOPE,
                                                    SLOPEDIST,
                                                    SHELFDIST,
                                                    SST,
                                                    SSTFRONT,
                                                    SSTVAR,
                                                    ICE,
                                                    ICEDIST,
                                                    ICEVAR)),
                    y = d2_Pacific$real,
                    method = "ranger",
                    metric = "ROC",
                    trControl = tc,
                    tuneGrid = rfGrid,
                    importance = "impurity",
                    num.trees = number_trees)

# Inspect
Mr_Pacific

# Varimp
varimpR(mod = Mr_Pacific, mod_name = "Mr_Pacific")

# Save
saveRDS(Mr_Pacific, "./output/fitted_models/Mr_Pacific.RDS")

# Performance
dat$Mr_Pacific_predicted_prob <- predict(Mr_Pacific, newdata = dat, type = "prob")$O
dat$Mr_Pacific_predicted_class <- predict(Mr_Pacific, newdata = dat, type = "raw")

internal_performance <- Mr_Pacific$results[which(Mr_Pacific$results$ROC == max(Mr_Pacific$results$ROC)),]$ROC
internal_performance_sd <- Mr_Pacific$results[which(Mr_Pacific$results$ROC == max(Mr_Pacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$Mr_Pacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$Mr_Pacific_p <- NA
area_grid$Mr_Pacific_p[dx] <- predict(Mr_Pacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/Mr_Pacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$Mr_Pacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "Mr_Pacific_p")])
writeRaster(r, "./output/predictions/Mr_Pacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- Mr_Pacific$results[which(Mr_Pacific$results$ROC == max(Mr_Pacific$results$ROC)),]
results$model <- "Mr_Pacific"
results$n_tracks <- length(unique(d2_Pacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_Mr_Pacific.csv", row.names = F)
rm(results)

#----------------------
# b. EastPacific

# Data
d2_EastPacific <- dplyr::filter(dat, region == "EastPacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR)

# Folds
folds <- groupKFold(group = d2_EastPacific$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
Mr_EastPacific <- train(x = as.data.frame(dplyr::select(d2_EastPacific,
                                                        DEPTH,
                                                        SLOPE,
                                                        SLOPEDIST,
                                                        SHELFDIST,
                                                        SST,
                                                        SSTFRONT,
                                                        SSTVAR,
                                                        ICE,
                                                        ICEDIST,
                                                        ICEVAR)),
                    y = d2_EastPacific$real,
                    method = "ranger",
                    metric = "ROC",
                    trControl = tc,
                    tuneGrid = rfGrid,
                    importance = "impurity",
                    num.trees = number_trees)

# Inspect
Mr_EastPacific

# Varimp
varimpR(mod = Mr_EastPacific, mod_name = "Mr_EastPacific")

# Save
saveRDS(Mr_EastPacific, "./output/fitted_models/Mr_EastPacific.RDS")

# Performance
dat$Mr_EastPacific_predicted_prob <- predict(Mr_EastPacific, newdata = dat, type = "prob")$O
dat$Mr_EastPacific_predicted_class <- predict(Mr_EastPacific, newdata = dat, type = "raw")

internal_performance <- Mr_EastPacific$results[which(Mr_EastPacific$results$ROC == max(Mr_EastPacific$results$ROC)),]$ROC
internal_performance_sd <- Mr_EastPacific$results[which(Mr_EastPacific$results$ROC == max(Mr_EastPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$Mr_EastPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$Mr_EastPacific_p <- NA
area_grid$Mr_EastPacific_p[dx] <- predict(Mr_EastPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/Mr_EastPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$Mr_EastPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "Mr_EastPacific_p")])
writeRaster(r, "./output/predictions/Mr_EastPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- Mr_EastPacific$results[which(Mr_EastPacific$results$ROC == max(Mr_EastPacific$results$ROC)),]
results$model <- "Mr_EastPacific"
results$n_tracks <- length(unique(d2_EastPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_Mr_EastPacific.csv", row.names = F)
rm(results)

#----------------------
# c. WestPacific

# Data
d2_WestPacific <- dplyr::filter(dat, region == "WestPacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR)

# Folds
folds <- groupKFold(group = d2_WestPacific$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
Mr_WestPacific <- train(x = as.data.frame(dplyr::select(d2_WestPacific,
                                                        DEPTH,
                                                        SLOPE,
                                                        SLOPEDIST,
                                                        SHELFDIST,
                                                        SST,
                                                        SSTFRONT,
                                                        SSTVAR,
                                                        ICE,
                                                        ICEDIST,
                                                        ICEVAR)),
                        y = d2_WestPacific$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                        importance = "impurity",
                        num.trees = number_trees)

# Inspect
Mr_WestPacific

# Varimp
varimpR(mod = Mr_WestPacific, mod_name = "Mr_WestPacific")

# Save
saveRDS(Mr_WestPacific, "./output/fitted_models/Mr_WestPacific.RDS")

# Performance
dat$Mr_WestPacific_predicted_prob <- predict(Mr_WestPacific, newdata = dat, type = "prob")$O
dat$Mr_WestPacific_predicted_class <- predict(Mr_WestPacific, newdata = dat, type = "raw")

internal_performance <- Mr_WestPacific$results[which(Mr_WestPacific$results$ROC == max(Mr_WestPacific$results$ROC)),]$ROC
internal_performance_sd <- Mr_WestPacific$results[which(Mr_WestPacific$results$ROC == max(Mr_WestPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$Mr_WestPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$Mr_WestPacific_p <- NA
area_grid$Mr_WestPacific_p[dx] <- predict(Mr_WestPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/Mr_WestPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$Mr_WestPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "Mr_WestPacific_p")])
writeRaster(r, "./output/predictions/Mr_WestPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- Mr_WestPacific$results[which(Mr_WestPacific$results$ROC == max(Mr_WestPacific$results$ROC)),]
results$model <- "Mr_WestPacific"
results$n_tracks <- length(unique(d2_WestPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_Mr_WestPacific.csv", row.names = F)
rm(results)

#----------------------
# d. EastIndian

# Data
d2_EastIndian <- dplyr::filter(dat, region == "EastIndian") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR)

# Folds
folds <- groupKFold(group = d2_EastIndian$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
Mr_EastIndian <- train(x = as.data.frame(dplyr::select(d2_EastIndian,
                                                       DEPTH,
                                                       SLOPE,
                                                       SLOPEDIST,
                                                       SHELFDIST,
                                                       SST,
                                                       SSTFRONT,
                                                       SSTVAR,
                                                       ICE,
                                                       ICEDIST,
                                                       ICEVAR)),
                        y = d2_EastIndian$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                       importance = "impurity",
                       num.trees = number_trees)

# Inspect
Mr_EastIndian

# Varimp
varimpR(mod = Mr_EastIndian, mod_name = "Mr_EastIndian")

# Save
saveRDS(Mr_EastIndian, "./output/fitted_models/Mr_EastIndian.RDS")

# Performance
dat$Mr_EastIndian_predicted_prob <- predict(Mr_EastIndian, newdata = dat, type = "prob")$O
dat$Mr_EastIndian_predicted_class <- predict(Mr_EastIndian, newdata = dat, type = "raw")

internal_performance <- Mr_EastIndian$results[which(Mr_EastIndian$results$ROC == max(Mr_EastIndian$results$ROC)),]$ROC
internal_performance_sd <- Mr_EastIndian$results[which(Mr_EastIndian$results$ROC == max(Mr_EastIndian$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$Mr_EastIndian_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$Mr_EastIndian_p <- NA
area_grid$Mr_EastIndian_p[dx] <- predict(Mr_EastIndian, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/Mr_EastIndian_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$Mr_EastIndian_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "Mr_EastIndian_p")])
writeRaster(r, "./output/predictions/Mr_EastIndian_raster.grd", overwrite = T)
rm(r)

# Save results
results <- Mr_EastIndian$results[which(Mr_EastIndian$results$ROC == max(Mr_EastIndian$results$ROC)),]
results$model <- "Mr_EastIndian"
results$n_tracks <- length(unique(d2_EastIndian$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_Mr_EastIndian.csv", row.names = F)
rm(results)

#----------------------
# e. Atlantic

# Data
d2_Atlantic <- dplyr::filter(dat, region == "Atlantic") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR)

# Folds
folds <- groupKFold(group = d2_Atlantic$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
Mr_Atlantic <- train(x = as.data.frame(dplyr::select(d2_Atlantic,
                                                     DEPTH,
                                                     SLOPE,
                                                     SLOPEDIST,
                                                     SHELFDIST,
                                                     SST,
                                                     SSTFRONT,
                                                     SSTVAR,
                                                     ICE,
                                                     ICEDIST,
                                                     ICEVAR)),
                        y = d2_Atlantic$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                        importance = "impurity")

# Inspect
Mr_Atlantic

# Varimp
varimpR(mod = Mr_Atlantic, mod_name = "Mr_Atlantic")

# Save
saveRDS(Mr_Atlantic, "./output/fitted_models/Mr_Atlantic.RDS")

# Performance
dat$Mr_Atlantic_predicted_prob <- predict(Mr_Atlantic, newdata = dat, type = "prob")$O
dat$Mr_Atlantic_predicted_class <- predict(Mr_Atlantic, newdata = dat, type = "raw")

internal_performance <- Mr_Atlantic$results[which(Mr_Atlantic$results$ROC == max(Mr_Atlantic$results$ROC)),]$ROC
internal_performance_sd <- Mr_Atlantic$results[which(Mr_Atlantic$results$ROC == max(Mr_Atlantic$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$Mr_Atlantic_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$Mr_Atlantic_p <- NA
area_grid$Mr_Atlantic_p[dx] <- predict(Mr_Atlantic, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/Mr_Atlantic_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$Mr_Atlantic_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "Mr_Atlantic_p")])
writeRaster(r, "./output/predictions/Mr_Atlantic_raster.grd", overwrite = T)
rm(r)

# Save results
results <- Mr_Atlantic$results[which(Mr_Atlantic$results$ROC == max(Mr_Atlantic$results$ROC)),]
results$model <- "Mr_Atlantic"
results$n_tracks <- length(unique(d2_Atlantic$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_Mr_Atlantic.csv", row.names = F)
rm(results)

#----------------------
# Create the ensemble prediction
dat$M2_predicted_prob <- rowMeans(dat[, c("Mr_Pacific_predicted_prob", "Mr_EastPacific_predicted_prob", "Mr_WestPacific_predicted_prob", "Mr_EastIndian_predicted_prob", "Mr_Atlantic_predicted_prob")])

# Track validation
performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_predicted_prob, levels = c("S", "O")))[1]

# IWC validation
area_grid$M2_p <- rowMeans(area_grid[, c("Mr_Pacific_p", "Mr_EastPacific_p", "Mr_WestPacific_p", "Mr_EastIndian_p", "Mr_Atlantic_p")])
saveRDS(area_grid, "./output/predictions/M2_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_p")])
writeRaster(r, "./output/predictions/M2_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- data.frame("mtry" = NA,
                      "splitrule" = NA,
                      "min.node.size" = NA,
                      "ROC" = NA,
                      "Sens" = NA,
                      "Spec" = NA,
                      "ROCSD" = NA,
                      "SensSD" = NA,
                      "SpecSD" = NA)
results$model <- "M2"
results$n_tracks <- length(unique(d1$id)) # Uses all the tracks
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2.csv", row.names = F)
rm(results)

#--------------------------------------------
# 3. Similarity-weighted ensemble
#--------------------------------------------
# Here we create a weighted ensemble of the regional predictions, training
# a model to predict membership of each datum to a region

# Data, use only real tracks
d3 <- dat[!is.na(dat$region), ]
d3 <- d3[d3$n_sim == 0, ]

folds <- groupKFold(group = d3$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   # sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = multiClassSummary,
                   index = folds)

# Fit the model
M3 <- train(x = as.data.frame(dplyr::select(d3,
                                            DEPTH,
                                            SLOPE,
                                            SLOPEDIST,
                                            SHELFDIST,
                                            SST,
                                            SSTFRONT,
                                            SSTVAR,
                                            ICE,
                                            ICEDIST,
                                            ICEVAR)),
            y = as.factor(d3$region),
            method = "ranger",
            metric = "AUC",
            trControl = tc,
            tuneGrid = rfGrid,
            importance = "impurity",
            num.trees = number_trees)

# Varimp
varimpR(mod = M3, mod_name = "M3")

# Save
saveRDS(M3, "./output/fitted_models/M3.RDS")

#----------------------
# Performance
internal_performance <- M3$results[which(M3$results$AUC == max(M3$results$AUC)),]$AUC

# Predict probability of a datum being in each region
dat$p_Atlantic <- predict(M3, newdata = dat, type = "prob")$Atlantic
dat$p_EastIndian <- predict(M3, newdata = dat, type = "prob")$EastIndian
dat$p_EastPacific <- predict(M3, newdata = dat, type = "prob")$EastPacific
dat$p_Pacific <- predict(M3, newdata = dat, type = "prob")$Pacific
dat$p_WestPacific <- predict(M3, newdata = dat, type = "prob")$WestPacific

# Multiply region membership probs by predictions from regional models
dat$M3_Atlantic_predicted_prob <- dat$Mr_Atlantic_predicted_prob * dat$p_Atlantic
dat$M3_EastIndian_predicted_prob <- dat$Mr_EastIndian_predicted_prob * dat$p_EastIndian
dat$M3_EastPacific_predicted_prob <- dat$Mr_EastPacific_predicted_prob * dat$p_EastPacific
dat$M3_Pacific_predicted_prob <- dat$Mr_Pacific_predicted_prob * dat$p_Pacific
dat$M3_WestPacific_predicted_prob <- dat$Mr_WestPacific_predicted_prob * dat$p_WestPacific

# Take mean
dat$M3_predicted_prob <- rowMeans(dat[, c("M3_Pacific_predicted_prob", "M3_EastPacific_predicted_prob", "M3_WestPacific_predicted_prob", "M3_EastIndian_predicted_prob", "M3_Atlantic_predicted_prob")])

# Track validation
performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M3_predicted_prob, levels = c("S", "O")))[1]


# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases

# Predict probability of a datum being in each region
area_grid$p_Atlantic <- NA
area_grid$p_EastIndian <- NA
area_grid$p_EastPacific <- NA
area_grid$p_Pacific <- NA
area_grid$p_WestPacific <- NA

area_grid$p_Atlantic[dx] <- predict(M3, newdata = area_grid[dx,], type = "prob")$Atlantic
area_grid$p_EastIndian[dx] <- predict(M3, newdata = area_grid[dx,], type = "prob")$EastIndian
area_grid$p_EastPacific[dx] <- predict(M3, newdata = area_grid[dx,], type = "prob")$EastPacific
area_grid$p_Pacific[dx] <- predict(M3, newdata = area_grid[dx,], type = "prob")$Pacific
area_grid$p_WestPacific[dx] <- predict(M3, newdata = area_grid[dx,], type = "prob")$WestPacific

# Multiply region membership probs by predictions from regional models
area_grid$M3_Atlantic_predicted_prob <- area_grid$Mr_Atlantic_p * area_grid$p_Atlantic
area_grid$M3_EastIndian_predicted_prob <- area_grid$Mr_EastIndian_p * area_grid$p_EastIndian
area_grid$M3_EastPacific_predicted_prob <- area_grid$Mr_EastPacific_p * area_grid$p_EastPacific
area_grid$M3_Pacific_predicted_prob <- area_grid$Mr_Pacific_p * area_grid$p_Pacific
area_grid$M3_WestPacific_predicted_prob <- area_grid$Mr_WestPacific_p * area_grid$p_WestPacific

# Take mean
area_grid$M3_predicted_prob <- rowMeans(area_grid[, c("M3_Pacific_predicted_prob", "M3_EastPacific_predicted_prob", "M3_WestPacific_predicted_prob", "M3_EastIndian_predicted_prob", "M3_Atlantic_predicted_prob")])

# Save
saveRDS(area_grid, "./output/predictions/M3_grid.RDS")

# IWC validation
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M3_predicted_prob, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M3_predicted_prob")])
writeRaster(r, "./output/predictions/M3_raster.grd", overwrite = T)
rm(r)


# Create a log file of results
results <- data.frame("mtry" = NA,
                      "splitrule" = NA,
                      "min.node.size" = NA,
                      "ROC" = NA,
                      "Sens" = NA,
                      "Spec" = NA,
                      "ROCSD" = NA,
                      "SensSD" = NA,
                      "SpecSD" = NA)
results$model <- "M3"
results$n_tracks <- length(unique(dat$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M3.csv", row.names = F)
rm(results)

#--------------------------------------------
# 4. Meta-model / simple ensemble - M4
#--------------------------------------------
# Here we use the predictions from the regional models to fit a new meta-model

# Data
d4 <- dat # Using everything

# Folds
folds <- groupKFold(group = d4$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)


# Fit the model
M4 <- train(x = as.data.frame(dplyr::select(d4,
                                            Mr_Pacific_predicted_prob,
                                            Mr_EastPacific_predicted_prob,
                                            Mr_WestPacific_predicted_prob,
                                            Mr_EastIndian_predicted_prob,
                                            Mr_Atlantic_predicted_prob)),
              y = d4$real,
              method = "ranger",
              metric = "ROC",
              trControl = tc,
              tuneGrid = rfGrid,
            importance = "impurity",
            num.trees = number_trees)

# Varimp
varimpR(mod = M4, mod_name = "M4")

# Save
saveRDS(M4, "./output/fitted_models/M4.RDS")

#-----------------
# Predict
# Performance
dat$M4_predicted_prob <- predict(M4, newdata = dat, type = "prob")$O
dat$M4_predicted_class <- predict(M4, newdata = dat, type = "raw")

internal_performance <- M4$results[which(M4$results$ROC == max(M4$results$ROC)),]$ROC
internal_performance_sd <- M4$results[which(M4$results$ROC == max(M4$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M4_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M4_p <- NA

# Create a temp grid with the same names
tmp <- mutate(area_grid,
              Mr_Pacific_predicted_prob = Mr_Pacific_p,
              Mr_EastPacific_predicted_prob = Mr_EastPacific_p,
              Mr_WestPacific_predicted_prob = Mr_WestPacific_p,
              Mr_EastIndian_predicted_prob = Mr_EastIndian_p,
              Mr_Atlantic_predicted_prob = Mr_Atlantic_p)
area_grid$M4_p[dx] <- predict(M4, newdata = tmp[dx, ], type = "prob")$O
rm(tmp)
saveRDS(area_grid, "./output/predictions/M4_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M4_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M4_p")])
writeRaster(r, "./output/predictions/M4_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M4$results[which(M4$results$ROC == max(M4$results$ROC)),]
results$model <- "M4"
results$n_tracks <- length(unique(d4$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M4.csv", row.names = F)
rm(results)

#--------------------------------------------
# 5. Informed ensemble - M5
#--------------------------------------------
# Here we include the predictions from the regional models to fit an informed model

# Data
d5 <- dat # Using everything

# Folds
folds <- groupKFold(group = d5$id, k = number_k)

# Train control
tc <- trainControl(method = "repeatedcv",
                   repeats = number_repeats,
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)


# Fit the model
M5 <- train(x = as.data.frame(dplyr::select(d5,
                                            DEPTH,
                                            SLOPE,
                                            SLOPEDIST,
                                            SHELFDIST,
                                            SST,
                                            SSTFRONT,
                                            SSTVAR,
                                            ICE,
                                            ICEDIST,
                                            ICEVAR,
                                            Mr_Pacific_predicted_prob,
                                            Mr_EastPacific_predicted_prob,
                                            Mr_WestPacific_predicted_prob,
                                            Mr_EastIndian_predicted_prob,
                                            Mr_Atlantic_predicted_prob)),
            y = d5$real,
            method = "ranger",
            metric = "ROC",
            trControl = tc,
            tuneGrid = rfGrid,
            importance = "impurity",
            num.trees = number_trees)

# Varimp
varimpR(mod = M5, mod_name = "M5")

# Save
saveRDS(M5, "./output/fitted_models/M5.RDS")

#-----------------
# Predict
# Performance
dat$M5_predicted_prob <- predict(M5, newdata = dat, type = "prob")$O
dat$M5_predicted_class <- predict(M5, newdata = dat, type = "raw")

internal_performance <- M5$results[which(M5$results$ROC == max(M5$results$ROC)),]$ROC
internal_performance_sd <- M5$results[which(M5$results$ROC == max(M5$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M5_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M5_p <- NA
# Create a temp grid with the same names
tmp <- mutate(area_grid,
              Mr_Pacific_predicted_prob = Mr_Pacific_p,
              Mr_EastPacific_predicted_prob = Mr_EastPacific_p,
              Mr_WestPacific_predicted_prob = Mr_WestPacific_p,
              Mr_EastIndian_predicted_prob = Mr_EastIndian_p,
              Mr_Atlantic_predicted_prob = Mr_Atlantic_p)
area_grid$M5_p[dx] <- predict(M5, newdata = tmp[dx, ], type = "prob")$O
rm(tmp)
saveRDS(area_grid, "./output/predictions/M5_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M5_p")])
writeRaster(r, "./output/predictions/M5_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M5$results[which(M5$results$ROC == max(M5$results$ROC)),]
results$model <- "M5"
results$n_tracks <- length(unique(d5$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M5.csv", row.names = F)
rm(results)