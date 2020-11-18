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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)

# Set up parameter grid
rfGrid <-  data.frame(mtry = c(3, 4, 5),
                      splitrule = "gini",
                      min.node.size = 1)

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
              importance = "impurity")
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
M2_Pacific <- train(x = as.data.frame(dplyr::select(d2_Pacific,
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
                    importance = "impurity")

# Inspect
M2_Pacific

# Varimp
varimpR(mod = M2_Pacific, mod_name = "M2_Pacific")

# Save
saveRDS(M2_Pacific, "./output/fitted_models/M2_Pacific.RDS")

# Performance
dat$M2_Pacific_predicted_prob <- predict(M2_Pacific, newdata = dat, type = "prob")$O
dat$M2_Pacific_predicted_class <- predict(M2_Pacific, newdata = dat, type = "raw")

internal_performance <- M2_Pacific$results[which(M2_Pacific$results$ROC == max(M2_Pacific$results$ROC)),]$ROC
internal_performance_sd <- M2_Pacific$results[which(M2_Pacific$results$ROC == max(M2_Pacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_Pacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_Pacific_p <- NA
area_grid$M2_Pacific_p[dx] <- predict(M2_Pacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M2_Pacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_Pacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_Pacific_p")])
writeRaster(r, "./output/predictions/M2_Pacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M2_Pacific$results[which(M2_Pacific$results$ROC == max(M2_Pacific$results$ROC)),]
results$model <- "M2_Pacific"
results$n_tracks <- length(unique(d2_Pacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2_Pacific.csv", row.names = F)
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
M2_EastPacific <- train(x = as.data.frame(dplyr::select(d2_EastPacific,
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
                    importance = "impurity")

# Inspect
M2_EastPacific

# Varimp
varimpR(mod = M2_EastPacific, mod_name = "M2_EastPacific")

# Save
saveRDS(M2_EastPacific, "./output/fitted_models/M2_EastPacific.RDS")

# Performance
dat$M2_EastPacific_predicted_prob <- predict(M2_EastPacific, newdata = dat, type = "prob")$O
dat$M2_EastPacific_predicted_class <- predict(M2_EastPacific, newdata = dat, type = "raw")

internal_performance <- M2_EastPacific$results[which(M2_EastPacific$results$ROC == max(M2_EastPacific$results$ROC)),]$ROC
internal_performance_sd <- M2_EastPacific$results[which(M2_EastPacific$results$ROC == max(M2_EastPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_EastPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_EastPacific_p <- NA
area_grid$M2_EastPacific_p[dx] <- predict(M2_EastPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M2_EastPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_EastPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_EastPacific_p")])
writeRaster(r, "./output/predictions/M2_EastPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M2_EastPacific$results[which(M2_EastPacific$results$ROC == max(M2_EastPacific$results$ROC)),]
results$model <- "M2_EastPacific"
results$n_tracks <- length(unique(d2_EastPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2_EastPacific.csv", row.names = F)
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
M2_WestPacific <- train(x = as.data.frame(dplyr::select(d2_WestPacific,
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
                        importance = "impurity")

# Inspect
M2_WestPacific

# Varimp
varimpR(mod = M2_WestPacific, mod_name = "M2_WestPacific")

# Save
saveRDS(M2_WestPacific, "./output/fitted_models/M2_WestPacific.RDS")

# Performance
dat$M2_WestPacific_predicted_prob <- predict(M2_WestPacific, newdata = dat, type = "prob")$O
dat$M2_WestPacific_predicted_class <- predict(M2_WestPacific, newdata = dat, type = "raw")

internal_performance <- M2_WestPacific$results[which(M2_WestPacific$results$ROC == max(M2_WestPacific$results$ROC)),]$ROC
internal_performance_sd <- M2_WestPacific$results[which(M2_WestPacific$results$ROC == max(M2_WestPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_WestPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_WestPacific_p <- NA
area_grid$M2_WestPacific_p[dx] <- predict(M2_WestPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M2_WestPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_WestPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_WestPacific_p")])
writeRaster(r, "./output/predictions/M2_WestPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M2_WestPacific$results[which(M2_WestPacific$results$ROC == max(M2_WestPacific$results$ROC)),]
results$model <- "M2_WestPacific"
results$n_tracks <- length(unique(d2_WestPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2_WestPacific.csv", row.names = F)
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
M2_EastIndian <- train(x = as.data.frame(dplyr::select(d2_EastIndian,
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
                        importance = "impurity")

# Inspect
M2_EastIndian

# Varimp
varimpR(mod = M2_EastIndian, mod_name = "M2_EastIndian")

# Save
saveRDS(M2_EastIndian, "./output/fitted_models/M2_EastIndian.RDS")

# Performance
dat$M2_EastIndian_predicted_prob <- predict(M2_EastIndian, newdata = dat, type = "prob")$O
dat$M2_EastIndian_predicted_class <- predict(M2_EastIndian, newdata = dat, type = "raw")

internal_performance <- M2_EastIndian$results[which(M2_EastIndian$results$ROC == max(M2_EastIndian$results$ROC)),]$ROC
internal_performance_sd <- M2_EastIndian$results[which(M2_EastIndian$results$ROC == max(M2_EastIndian$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_EastIndian_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_EastIndian_p <- NA
area_grid$M2_EastIndian_p[dx] <- predict(M2_EastIndian, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M2_EastIndian_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_EastIndian_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_EastIndian_p")])
writeRaster(r, "./output/predictions/M2_EastIndian_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M2_EastIndian$results[which(M2_EastIndian$results$ROC == max(M2_EastIndian$results$ROC)),]
results$model <- "M2_EastIndian"
results$n_tracks <- length(unique(d2_EastIndian$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2_EastIndian.csv", row.names = F)
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)
# Fit
M2_Atlantic <- train(x = as.data.frame(dplyr::select(d2_Atlantic,
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
M2_Atlantic

# Varimp
varimpR(mod = M2_Atlantic, mod_name = "M2_Atlantic")

# Save
saveRDS(M2_Atlantic, "./output/fitted_models/M2_Atlantic.RDS")

# Performance
dat$M2_Atlantic_predicted_prob <- predict(M2_Atlantic, newdata = dat, type = "prob")$O
dat$M2_Atlantic_predicted_class <- predict(M2_Atlantic, newdata = dat, type = "raw")

internal_performance <- M2_Atlantic$results[which(M2_Atlantic$results$ROC == max(M2_Atlantic$results$ROC)),]$ROC
internal_performance_sd <- M2_Atlantic$results[which(M2_Atlantic$results$ROC == max(M2_Atlantic$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_Atlantic_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_Atlantic_p <- NA
area_grid$M2_Atlantic_p[dx] <- predict(M2_Atlantic, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M2_Atlantic_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_Atlantic_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_Atlantic_p")])
writeRaster(r, "./output/predictions/M2_Atlantic_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M2_Atlantic$results[which(M2_Atlantic$results$ROC == max(M2_Atlantic$results$ROC)),]
results$model <- "M2_Atlantic"
results$n_tracks <- length(unique(d2_Atlantic$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2_Atlantic.csv", row.names = F)
rm(results)

#----------------------
# Create the ensemble prediction
dat$M2_predicted_prob <- rowMeans(dat[, c("M2_Pacific_predicted_prob", "M2_EastPacific_predicted_prob", "M2_WestPacific_predicted_prob", "M2_EastIndian_predicted_prob", "M2_Atlantic_predicted_prob")])

# Track validation
performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_predicted_prob, levels = c("S", "O")))[1]

# IWC validation
area_grid$M2_p <- rowMeans(area_grid[, c("M2_Pacific_p", "M2_EastPacific_p", "M2_WestPacific_p", "M2_EastIndian_p", "M2_Atlantic_p")])
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
results$n_tracks <- length(unique(d3$id))
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
tc <- trainControl(method = "cv",
                   number = number_k,
                   search = "grid",
                   classProbs = TRUE,
                   # sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = multiClassSummary,
                   index = folds)

# Params
rfGrid <-  data.frame(mtry = c(3, 4, 5),
                      splitrule = "gini",
                      min.node.size = 1)

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
            importance = "impurity")

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
dat$M3_Atlantic_predicted_prob <- dat$M2_Atlantic_predicted_prob * dat$p_Atlantic
dat$M3_EastIndian_predicted_prob <- dat$M2_EastIndian_predicted_prob * dat$p_EastIndian
dat$M3_EastPacific_predicted_prob <- dat$M2_EastPacific_predicted_prob * dat$p_EastPacific
dat$M3_Pacific_predicted_prob <- dat$M2_Pacific_predicted_prob * dat$p_Pacific
dat$M3_WestPacific_predicted_prob <- dat$M2_WestPacific_predicted_prob * dat$p_WestPacific

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
area_grid$M3_Atlantic_predicted_prob <- area_grid$M2_Atlantic_p * area_grid$p_Atlantic
area_grid$M3_EastIndian_predicted_prob <- area_grid$M2_EastIndian_p * area_grid$p_EastIndian
area_grid$M3_EastPacific_predicted_prob <- area_grid$M2_EastPacific_p * area_grid$p_EastPacific
area_grid$M3_Pacific_predicted_prob <- area_grid$M2_Pacific_p * area_grid$p_Pacific
area_grid$M3_WestPacific_predicted_prob <- area_grid$M2_WestPacific_p * area_grid$p_WestPacific

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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)

# Params
rfGrid <-  data.frame(mtry = c(3, 4, 5),
                      splitrule = "gini",
                      min.node.size = 1)

# Fit the model
M4 <- train(x = as.data.frame(dplyr::select(d4,
                                            M2_Pacific_predicted_prob,
                                            M2_EastPacific_predicted_prob,
                                            M2_WestPacific_predicted_prob,
                                            M2_EastIndian_predicted_prob,
                                            M2_Atlantic_predicted_prob)),
              y = d4$real,
              method = "ranger",
              metric = "ROC",
              trControl = tc,
              tuneGrid = rfGrid,
              importance = "impurity")

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
              M2_Pacific_predicted_prob = M2_Pacific_p,
              M2_EastPacific_predicted_prob = M2_EastPacific_p,
              M2_WestPacific_predicted_prob = M2_WestPacific_p,
              M2_EastIndian_predicted_prob = M2_EastIndian_p,
              M2_Atlantic_predicted_prob = M2_Atlantic_p)
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
tc <- trainControl(method = "cv",
                   number = length(folds),
                   search = "grid",
                   classProbs = TRUE,
                   sampling = "down",
                   allowParallel = allow_par,
                   summaryFunction = twoClassSummary,
                   index = folds)

# Params
rfGrid <-  data.frame(mtry = c(3, 4, 5),
                      splitrule = "gini",
                      min.node.size = 1)

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
                                            M2_Pacific_predicted_prob,
                                            M2_EastPacific_predicted_prob,
                                            M2_WestPacific_predicted_prob,
                                            M2_EastIndian_predicted_prob,
                                            M2_Atlantic_predicted_prob)),
            y = d5$real,
            method = "ranger",
            metric = "ROC",
            trControl = tc,
            tuneGrid = rfGrid,
            importance = "impurity")

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
              M2_Pacific_predicted_prob = M2_Pacific_p,
              M2_EastPacific_predicted_prob = M2_EastPacific_p,
              M2_WestPacific_predicted_prob = M2_WestPacific_p,
              M2_EastIndian_predicted_prob = M2_EastIndian_p,
              M2_Atlantic_predicted_prob = M2_Atlantic_p)
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