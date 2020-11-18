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

# Depending on the model, caret can only run with complete cases
dat <- dat[complete.cases(dat[,9:21]), ]

#-------------------------------------------------------------------
# Load data for prediction
area_grid <- readRDS("./output/envar_frame.RDS")

# Add IWC data for validation
iwc <- raster("./output/iwc_catches_and_survey.grd")
area_grid$iwc <- raster::extract(x = iwc, y = area_grid[,c("x", "y")])
area_grid$iwc[which(area_grid$iwc == 1)] <- "O"
area_grid$iwc[which(area_grid$iwc == 0)] <- "S"

#-------------------------------------------------------------------
# Function to create folds for individual based CV

foldCaret <- function(dat, nm = 10) {
  
  #--------------------------------------------------------
  # Create a folds index to keep cells together during cross-validation
  # The resulting list "folds" is passed to "index" in the "trainControl" function in caret
  # Keeps individuals together during splitting data for cross validation
  
  use.long <- dat
  nm <- nm
  
  if (length(unique(use.long$track_id)) < nm) {
    nm <- length(unique(use.long$track_id))
  }
  
  rws <- floor(length(unique(use.long$track_id))/nm) #number of cells to include in each fold
  cells <- unique(use.long$track_id)
  folds <- list()
  for (i in 1:nm) {
    samp <- sample(x = cells, size = rws)
    l <- list()
    for (j in 1:length(samp)) {
      k <- which(use.long$track_id == samp[j])
      l[j] <- list(k)
    }
    l <- unlist(l)
    folds[i] <- list(l)
    for (r in 1:length(samp)) {
      cells <- cells[cells != samp[r]]
    }
  }
  
  names(folds) <- paste0("fold_", c(1:nm)) #train needs list to be named
  
  return(folds)
  
}


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
# 1. The naive global model - M0
# Here we use all data, fit a single model, and predict to the whole area

# Select our data
d0 <- dat # Using everything

# Create folds
# folds <- foldCaret(dat = d0, nm = 10)
folds <- groupKFold(group = d0$id, k = number_k) # This is caret's new built-in function

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
  M0 <- train(x = as.data.frame(d0[ ,c(9:21)]),
              y = d0$real,
              method = "ranger",
              metric = "ROC",
              trControl = tc,
              tuneGrid = rfGrid,
              importance = "impurity")
)

# Inspect
M0

# Variable importance
varImp(M0)

# Run function to plot and save varimp
varimpR(mod = M0, mod_name = "M0")

# Save
saveRDS(M0, "./output/fitted_models/M0.RDS")

# Stop the cluster
# stopCluster(clust)
# registerDoSEQ()

#-----------------
# Predict
# Performance
dat$M0_predicted_prob <- predict(M0, newdata = dat, type = "prob")$O
dat$M0_predicted_class <- predict(M0, newdata = dat, type = "raw")

internal_performance <- M0$results[which(M0$results$ROC == max(M0$results$ROC)),]$ROC
internal_performance_sd <- M0$results[which(M0$results$ROC == max(M0$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M0_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M0_p <- NA
area_grid$M0_p[dx] <- predict(M0, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M0_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M0_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M0_p")])
writeRaster(r, "./output/predictions/M0_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M0$results[which(M0$results$ROC == max(M0$results$ROC)),]
results$model <- "M0"
results$n_tracks <- length(unique(d0$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M0.csv", row.names = F)
rm(results)
#--------------------------------------------

#--------------------------------------------
# 2. Unweighted ensemble / simple mean - M1
# Here we fit models to each region, and for the global area take the mean of their predictions

# First, we fit the regional models
#----------------------
# a. Pacific

# Data
d1_Pacific <- dplyr::filter(dat, region == "Pacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
         DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR, EKE, SSHA, SSHGRAD)

# Folds
folds <- groupKFold(group = d1_Pacific$id, k = number_k)

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
M1_Pacific <- train(x = as.data.frame(d1_Pacific[ ,c(9:21)]),
                    y = d1_Pacific$real,
                    method = "ranger",
                    metric = "ROC",
                    trControl = tc,
                    tuneGrid = rfGrid,
                    importance = "impurity")

# Inspect
M1_Pacific

# Varimp
varimpR(mod = M1_Pacific, mod_name = "M1_Pacific")

# Save
saveRDS(M1_Pacific, "./output/fitted_models/M1_Pacific.RDS")

# Performance
dat$M1_Pacific_predicted_prob <- predict(M1_Pacific, newdata = dat, type = "prob")$O
dat$M1_Pacific_predicted_class <- predict(M1_Pacific, newdata = dat, type = "raw")

internal_performance <- M1_Pacific$results[which(M1_Pacific$results$ROC == max(M1_Pacific$results$ROC)),]$ROC
internal_performance_sd <- M1_Pacific$results[which(M1_Pacific$results$ROC == max(M1_Pacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_Pacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_Pacific_p <- NA
area_grid$M1_Pacific_p[dx] <- predict(M1_Pacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_Pacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_Pacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_Pacific_p")])
writeRaster(r, "./output/predictions/M1_Pacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M1_Pacific$results[which(M1_Pacific$results$ROC == max(M1_Pacific$results$ROC)),]
results$model <- "M1_Pacific"
results$n_tracks <- length(unique(d1_Pacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1_Pacific.csv", row.names = F)
rm(results)

#----------------------
# b. EastPacific

# Data
d1_EastPacific <- dplyr::filter(dat, region == "EastPacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR, EKE, SSHA, SSHGRAD)

# Folds
folds <- groupKFold(group = d1_EastPacific$id, k = number_k)

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
M1_EastPacific <- train(x = as.data.frame(d1_EastPacific[ ,c(9:21)]),
                    y = d1_EastPacific$real,
                    method = "ranger",
                    metric = "ROC",
                    trControl = tc,
                    tuneGrid = rfGrid,
                    importance = "impurity")

# Inspect
M1_EastPacific

# Varimp
varimpR(mod = M1_EastPacific, mod_name = "M1_EastPacific")

# Save
saveRDS(M1_EastPacific, "./output/fitted_models/M1_EastPacific.RDS")

# Performance
dat$M1_EastPacific_predicted_prob <- predict(M1_EastPacific, newdata = dat, type = "prob")$O
dat$M1_EastPacific_predicted_class <- predict(M1_EastPacific, newdata = dat, type = "raw")

internal_performance <- M1_EastPacific$results[which(M1_EastPacific$results$ROC == max(M1_EastPacific$results$ROC)),]$ROC
internal_performance_sd <- M1_EastPacific$results[which(M1_EastPacific$results$ROC == max(M1_EastPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_EastPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_EastPacific_p <- NA
area_grid$M1_EastPacific_p[dx] <- predict(M1_EastPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_EastPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_EastPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_EastPacific_p")])
writeRaster(r, "./output/predictions/M1_EastPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M1_EastPacific$results[which(M1_EastPacific$results$ROC == max(M1_EastPacific$results$ROC)),]
results$model <- "M1_EastPacific"
results$n_tracks <- length(unique(d1_EastPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1_EastPacific.csv", row.names = F)
rm(results)

#----------------------
# c. WestPacific

# Data
d1_WestPacific <- dplyr::filter(dat, region == "WestPacific") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR, EKE, SSHA, SSHGRAD)

# Folds
folds <- groupKFold(group = d1_WestPacific$id, k = number_k)

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
M1_WestPacific <- train(x = as.data.frame(d1_WestPacific[ ,c(9:21)]),
                        y = d1_WestPacific$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                        importance = "impurity")

# Inspect
M1_WestPacific

# Varimp
varimpR(mod = M1_WestPacific, mod_name = "M1_WestPacific")

# Save
saveRDS(M1_WestPacific, "./output/fitted_models/M1_WestPacific.RDS")

# Performance
dat$M1_WestPacific_predicted_prob <- predict(M1_WestPacific, newdata = dat, type = "prob")$O
dat$M1_WestPacific_predicted_class <- predict(M1_WestPacific, newdata = dat, type = "raw")

internal_performance <- M1_WestPacific$results[which(M1_WestPacific$results$ROC == max(M1_WestPacific$results$ROC)),]$ROC
internal_performance_sd <- M1_WestPacific$results[which(M1_WestPacific$results$ROC == max(M1_WestPacific$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_WestPacific_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_WestPacific_p <- NA
area_grid$M1_WestPacific_p[dx] <- predict(M1_WestPacific, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_WestPacific_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_WestPacific_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_WestPacific_p")])
writeRaster(r, "./output/predictions/M1_WestPacific_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M1_WestPacific$results[which(M1_WestPacific$results$ROC == max(M1_WestPacific$results$ROC)),]
results$model <- "M1_WestPacific"
results$n_tracks <- length(unique(d1_WestPacific$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1_WestPacific.csv", row.names = F)
rm(results)

#----------------------
# d. EastIndian

# Data
d1_EastIndian <- dplyr::filter(dat, region == "EastIndian") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR, EKE, SSHA, SSHGRAD)

# Folds
folds <- groupKFold(group = d1_EastIndian$id, k = number_k)

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
M1_EastIndian <- train(x = as.data.frame(d1_EastIndian[ ,c(9:21)]),
                        y = d1_EastIndian$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                        importance = "impurity")

# Inspect
M1_EastIndian

# Varimp
varimpR(mod = M1_EastIndian, mod_name = "M1_EastIndian")

# Save
saveRDS(M1_EastIndian, "./output/fitted_models/M1_EastIndian.RDS")

# Performance
dat$M1_EastIndian_predicted_prob <- predict(M1_EastIndian, newdata = dat, type = "prob")$O
dat$M1_EastIndian_predicted_class <- predict(M1_EastIndian, newdata = dat, type = "raw")

internal_performance <- M1_EastIndian$results[which(M1_EastIndian$results$ROC == max(M1_EastIndian$results$ROC)),]$ROC
internal_performance_sd <- M1_EastIndian$results[which(M1_EastIndian$results$ROC == max(M1_EastIndian$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_EastIndian_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_EastIndian_p <- NA
area_grid$M1_EastIndian_p[dx] <- predict(M1_EastIndian, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_EastIndian_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_EastIndian_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_EastIndian_p")])
writeRaster(r, "./output/predictions/M1_EastIndian_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M1_EastIndian$results[which(M1_EastIndian$results$ROC == max(M1_EastIndian$results$ROC)),]
results$model <- "M1_EastIndian"
results$n_tracks <- length(unique(d1_EastIndian$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1_EastIndian.csv", row.names = F)
rm(results)

#----------------------
# e. Atlantic

# Data
d1_Atlantic <- dplyr::filter(dat, region == "Atlantic") %>% 
  dplyr::select(., id, date, lon, lat, n_sim, real, track_id, region,
                DEPTH, SLOPE, SLOPEDIST, SHELFDIST, SST, SSTFRONT, SSTVAR, ICE, ICEDIST, ICEVAR, EKE, SSHA, SSHGRAD)

# Folds
folds <- groupKFold(group = d1_Atlantic$id, k = number_k)

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
M1_Atlantic <- train(x = as.data.frame(d1_Atlantic[ ,c(9:21)]),
                        y = d1_Atlantic$real,
                        method = "ranger",
                        metric = "ROC",
                        trControl = tc,
                        tuneGrid = rfGrid,
                        importance = "impurity")

# Inspect
M1_Atlantic

# Varimp
varimpR(mod = M1_Atlantic, mod_name = "M1_Atlantic")

# Save
saveRDS(M1_Atlantic, "./output/fitted_models/M1_Atlantic.RDS")

# Performance
dat$M1_Atlantic_predicted_prob <- predict(M1_Atlantic, newdata = dat, type = "prob")$O
dat$M1_Atlantic_predicted_class <- predict(M1_Atlantic, newdata = dat, type = "raw")

internal_performance <- M1_Atlantic$results[which(M1_Atlantic$results$ROC == max(M1_Atlantic$results$ROC)),]$ROC
internal_performance_sd <- M1_Atlantic$results[which(M1_Atlantic$results$ROC == max(M1_Atlantic$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_Atlantic_predicted_prob, levels = c("S", "O")))[1]

dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M1_Atlantic_p <- NA
area_grid$M1_Atlantic_p[dx] <- predict(M1_Atlantic, newdata = area_grid[dx, ], type = "prob")$O
saveRDS(area_grid, "./output/predictions/M1_Atlantic_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_Atlantic_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_Atlantic_p")])
writeRaster(r, "./output/predictions/M1_Atlantic_raster.grd", overwrite = T)
rm(r)

# Save results
results <- M1_Atlantic$results[which(M1_Atlantic$results$ROC == max(M1_Atlantic$results$ROC)),]
results$model <- "M1_Atlantic"
results$n_tracks <- length(unique(d1_Atlantic$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M1_Atlantic.csv", row.names = F)
rm(results)

#----------------------
# Create the ensemble prediction
dat$M1_predicted_prob <- rowMeans(dat[, c("M1_Pacific_predicted_prob", "M1_EastPacific_predicted_prob", "M1_WestPacific_predicted_prob", "M1_EastIndian_predicted_prob", "M1_Atlantic_predicted_prob")])

# Track validation
performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M1_predicted_prob, levels = c("S", "O")))[1]

# IWC validation
area_grid$M1_p <- rowMeans(area_grid[, c("M1_Pacific_p", "M1_EastPacific_p", "M1_WestPacific_p", "M1_EastIndian_p", "M1_Atlantic_p")])
saveRDS(area_grid, "./output/predictions/M1_grid.RDS")

validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M1_p, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M1_p")])
writeRaster(r, "./output/predictions/M1_raster.grd", overwrite = T)
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
results$model <- "M1"
results$n_tracks <- length(unique(d3$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M3.csv", row.names = F)
rm(results)

#--------------------------------------------
# 3. Meta-model / simple ensemble - M2
# Here we use the predictions from the regional models to fit a new meta-model

# Data
d2 <- dat # Using everything

# Folds
folds <- groupKFold(group = d2$id, k = number_k)

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
M2 <- train(x = as.data.frame(dplyr::select(d2,
                                            M1_Pacific_predicted_prob,
                                            M1_EastPacific_predicted_prob,
                                            M1_WestPacific_predicted_prob,
                                            M1_EastIndian_predicted_prob,
                                            M1_Atlantic_predicted_prob)),
              y = d2$real,
              method = "ranger",
              metric = "ROC",
              trControl = tc,
              tuneGrid = rfGrid,
              importance = "impurity")

# Varimp
varimpR(mod = M2, mod_name = "M2")

# Save
saveRDS(M2, "./output/fitted_models/M2.RDS")

#-----------------
# Predict
# Performance
dat$M2_predicted_prob <- predict(M2, newdata = dat, type = "prob")$O
dat$M2_predicted_class <- predict(M2, newdata = dat, type = "raw")

internal_performance <- M2$results[which(M2$results$ROC == max(M2$results$ROC)),]$ROC
internal_performance_sd <- M2$results[which(M2$results$ROC == max(M2$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M2_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M2_p <- NA
# Create a temp grid with the same names
tmp <- mutate(area_grid,
              M1_Pacific_predicted_prob = M1_Pacific_p,
              M1_EastPacific_predicted_prob = M1_EastPacific_p,
              M1_WestPacific_predicted_prob = M1_WestPacific_p,
              M1_EastIndian_predicted_prob = M1_EastIndian_p,
              M1_Atlantic_predicted_prob = M1_Atlantic_p)
area_grid$M2_p[dx] <- predict(M2, newdata = tmp[dx, ], type = "prob")$O
rm(tmp)
saveRDS(area_grid, "./output/predictions/M2_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M2_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_p")])
writeRaster(r, "./output/predictions/M2_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M2$results[which(M2$results$ROC == max(M2$results$ROC)),]
results$model <- "M2"
results$n_tracks <- length(unique(d2$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M2.csv", row.names = F)
rm(results)

#--------------------------------------------
# 4. Informed ensemble
# Here we include the predictions from the regional models to fit an informed model

# Data
d3 <- dat # Using everything

# Folds
folds <- groupKFold(group = d3$id, k = number_k)

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
                                            ICEVAR,
                                            EKE,
                                            SSHA,
                                            SSHGRAD,
                                            M1_Pacific_predicted_prob,
                                            M1_EastPacific_predicted_prob,
                                            M1_WestPacific_predicted_prob,
                                            M1_EastIndian_predicted_prob,
                                            M1_Atlantic_predicted_prob)),
            y = d3$real,
            method = "ranger",
            metric = "ROC",
            trControl = tc,
            tuneGrid = rfGrid,
            importance = "impurity")

# Varimp
varimpR(mod = M3, mod_name = "M3")

# Save
saveRDS(M3, "./output/fitted_models/M3.RDS")

#-----------------
# Predict
# Performance
dat$M3_predicted_prob <- predict(M3, newdata = dat, type = "prob")$O
dat$M3_predicted_class <- predict(M3, newdata = dat, type = "raw")

internal_performance <- M3$results[which(M3$results$ROC == max(M3$results$ROC)),]$ROC
internal_performance_sd <- M3$results[which(M3$results$ROC == max(M3$results$ROC)),]$ROCSD

performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M3_predicted_prob, levels = c("S", "O")))[1]

# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases
area_grid$M3_p <- NA
# Create a temp grid with the same names
tmp <- mutate(area_grid,
              M1_Pacific_predicted_prob = M1_Pacific_p,
              M1_EastPacific_predicted_prob = M1_EastPacific_p,
              M1_WestPacific_predicted_prob = M1_WestPacific_p,
              M1_EastIndian_predicted_prob = M1_EastIndian_p,
              M1_Atlantic_predicted_prob = M1_Atlantic_p)
area_grid$M3_p[dx] <- predict(M3, newdata = tmp[dx, ], type = "prob")$O
rm(tmp)
saveRDS(area_grid, "./output/predictions/M3_grid.RDS")

# Validation with IWC
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M3_p, levels = c("S", "O")))[1]

#Create raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M2_p")])
writeRaster(r, "./output/predictions/M3_raster.grd", overwrite = T)
rm(r)

# Create a log file of results
results <- M3$results[which(M3$results$ROC == max(M3$results$ROC)),]
results$model <- "M3"
results$n_tracks <- length(unique(d3$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M3.csv", row.names = F)
rm(results)

#--------------------------------------------
# 5. Weighted ensemble
# Here we create a weighted ensemble of the regional predictions, training
# a model to predict membership of each datum to a region

# Data, use only real tracks
d4 <- dat[!is.na(dat$region), ]
d4 <- d4[d4$n_sim == 0, ]

folds <- groupKFold(group = d4$id, k = number_k)

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
M4 <- train(x = as.data.frame(dplyr::select(d4,
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
                                            EKE,
                                            SSHA,
                                            SSHGRAD)),
            y = as.factor(d4$region),
            method = "ranger",
            metric = "AUC",
            trControl = tc,
            tuneGrid = rfGrid,
            importance = "impurity")

# Varimp
varimpR(mod = M4, mod_name = "M4")

# Save
saveRDS(M4, "./output/fitted_models/M4.RDS")

#----------------------
# Performance
internal_performance <- M4$results[which(M4$results$AUC == max(M4$results$AUC)),]$AUC

# Predict probability of a datum being in each region
dat$p_Atlantic <- predict(M4, newdata = dat, type = "prob")$Atlantic
dat$p_EastIndian <- predict(M4, newdata = dat, type = "prob")$EastIndian
dat$p_EastPacific <- predict(M4, newdata = dat, type = "prob")$EastPacific
dat$p_Pacific <- predict(M4, newdata = dat, type = "prob")$Pacific
dat$p_WestPacific <- predict(M4, newdata = dat, type = "prob")$WestPacific

# Multiply region membership probs by predictions from regional models
dat$M4_Atlantic_predicted_prob <- dat$M1_Atlantic_predicted_prob * dat$p_Atlantic
dat$M4_EastIndian_predicted_prob <- dat$M1_EastIndian_predicted_prob * dat$p_EastIndian
dat$M4_EastPacific_predicted_prob <- dat$M1_EastPacific_predicted_prob * dat$p_EastPacific
dat$M4_Pacific_predicted_prob <- dat$M1_Pacific_predicted_prob * dat$p_Pacific
dat$M4_WestPacific_predicted_prob <- dat$M1_WestPacific_predicted_prob * dat$p_WestPacific

# Take mean
dat$M4_predicted_prob <- rowMeans(dat[, c("M4_Pacific_predicted_prob", "M4_EastPacific_predicted_prob", "M4_WestPacific_predicted_prob", "M4_EastIndian_predicted_prob", "M4_Atlantic_predicted_prob")])

# Track validation
performance <- pROC::auc(pROC::roc(response = dat$real, predictor = dat$M4_predicted_prob, levels = c("S", "O")))[1]


# Predict to grid
dx <- which(complete.cases(area_grid)) # indicator for complete cases

# Predict probability of a datum being in each region
area_grid$p_Atlantic <- NA
area_grid$p_EastIndian <- NA
area_grid$p_EastPacific <- NA
area_grid$p_Pacific <- NA
area_grid$p_WestPacific <- NA

area_grid$p_Atlantic[dx] <- predict(M4, newdata = area_grid[dx,], type = "prob")$Atlantic
area_grid$p_EastIndian[dx] <- predict(M4, newdata = area_grid[dx,], type = "prob")$EastIndian
area_grid$p_EastPacific[dx] <- predict(M4, newdata = area_grid[dx,], type = "prob")$EastPacific
area_grid$p_Pacific[dx] <- predict(M4, newdata = area_grid[dx,], type = "prob")$Pacific
area_grid$p_WestPacific[dx] <- predict(M4, newdata = area_grid[dx,], type = "prob")$WestPacific

# Multiply region membership probs by predictions from regional models
area_grid$M4_Atlantic_predicted_prob <- area_grid$M1_Atlantic_p * area_grid$p_Atlantic
area_grid$M4_EastIndian_predicted_prob <- area_grid$M1_EastIndian_p * area_grid$p_EastIndian
area_grid$M4_EastPacific_predicted_prob <- area_grid$M1_EastPacific_p * area_grid$p_EastPacific
area_grid$M4_Pacific_predicted_prob <- area_grid$M1_Pacific_p * area_grid$p_Pacific
area_grid$M4_WestPacific_predicted_prob <- area_grid$M1_WestPacific_p * area_grid$p_WestPacific

# Take mean
area_grid$M4_predicted_prob <- rowMeans(area_grid[, c("M4_Pacific_predicted_prob", "M4_EastPacific_predicted_prob", "M4_WestPacific_predicted_prob", "M4_EastIndian_predicted_prob", "M4_Atlantic_predicted_prob")])

# Save
saveRDS(area_grid, "./output/predictions/M4_grid.RDS")

# IWC validation
validation_performance <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M4_predicted_prob, levels = c("S", "O")))[1]

# Save raster
r <- rasterFromXYZ(area_grid[,c("x", "y", "M4_predicted_prob")])
writeRaster(r, "./output/predictions/M4_raster.grd", overwrite = T)
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
results$model <- "M4"
results$n_tracks <- length(unique(dat$id))
results$validation_all_tracks <- performance
results$validation_iwc <- validation_performance
write.csv(results, "./output/fitted_models_results/results_M4.csv", row.names = F)
rm(results)