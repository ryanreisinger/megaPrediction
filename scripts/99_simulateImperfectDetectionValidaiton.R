setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

library(pROC)

area_grid <- readRDS("./output/predictions/M5_grid.RDS")

this_roc <- pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O"))
plot(this_roc)
this_auc <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O")))[1]

#---------------------------------------
# Lower the detection rate
plot(this_roc)
these_aucs_lower <- data.frame("iter" = 1:100,
                         "auc" = NA)
for (i in 1:50) {
  print(i)
area_grid$iwc_lower_detect <- area_grid$iwc
dx <- nrow(area_grid[area_grid$iwc_lower_detect == "O", ])
# dx <- sample(dx, size = floor(0.2*dx))
dx <- sample(dx, size = 50000)
area_grid[area_grid$iwc_lower_detect == "O", ]$iwc_lower_detect[dx] <- "S"

# Calculate 
this_roc_lower <- pROC::roc(response = area_grid$iwc_lower_detect, predictor = area_grid$M5_p, levels = c("S", "O"))
plot(this_roc_lower, add = T, col = "grey")
this_auc_lower <- pROC::auc(pROC::roc(response = area_grid$iwc_lower_detect, predictor = area_grid$M5_p, levels = c("S", "O")))[1]

# Add to output
these_aucs_lower$auc[i] <- this_auc_lower
}

#---------------------------------------
# Increase the detection rate
plot(this_roc)
these_aucs_higher <- data.frame("iter" = 1:100,
                               "auc" = NA)
for (i in 1:50) {
  print(i)
  area_grid$iwc_higher_detect <- area_grid$iwc
  dx <- nrow(area_grid[area_grid$iwc_higher_detect == "S", ])
  # dx <- sample(dx, size = floor(0.2*dx))
  dx <- sample(dx, size = 50000)
  area_grid[area_grid$iwc_higher_detect == "S", ]$iwc_higher_detect[dx] <- "O"
  
  # Calculate 
  this_roc_higher <- pROC::roc(response = area_grid$iwc_higher_detect, predictor = area_grid$M5_p, levels = c("S", "O"))
  plot(this_roc_higher, add = T, col = "grey")
  this_auc_higher <- pROC::auc(pROC::roc(response = area_grid$iwc_higher_detect, predictor = area_grid$M5_p, levels = c("S", "O")))[1]
  
  # Add to output
  these_aucs_higher$auc[i] <- this_auc_higher
}

#---------------------------------------
# Simulate true positives becoming false positives
plot(this_roc)
these_aucs_TP_FP <- data.frame("iter" = 1:10,
                                "auc" = NA)

d <- which(!is.na(area_grid$iwc) & !is.na(area_grid$M5_p) & area_grid$iwc == "O" & area_grid$M5_p > 0.75)

for (i in 1:10) {
  print(i)
  dx <- length(d)
  dx <- sample(d, size = floor(0.2*dx))
  area_grid$iwc[dx] <- "S"
  
  # Calculate 
  this_roc_higher <- pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O"))
  plot(this_roc_higher, add = T, col = "grey")
  this_auc_higher <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O")))[1]
  
  # Add to output
  these_aucs_TP_FP$auc[i] <- this_auc_higher
}

#---------------------------------------
# Simulate false negatives becoming true negatives
plot(this_roc)
these_aucs_FN_TN <- data.frame("iter" = 1:10,
                               "auc" = NA)

d <- which(!is.na(area_grid$iwc) & !is.na(area_grid$M5_p) & area_grid$iwc == "S" & area_grid$M5_p < 0.25)

for (i in 1:10) {
  print(i)
  dx <- length(d)
  dx <- sample(d, size = floor(0.2*dx))
  area_grid$iwc[dx] <- "O"
  
  # Calculate 
  this_roc_higher <- pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O"))
  plot(this_roc_higher, add = T, col = "grey")
  this_auc_higher <- pROC::auc(pROC::roc(response = area_grid$iwc, predictor = area_grid$M5_p, levels = c("S", "O")))[1]
  
  # Add to output
  these_aucs_FN_TN$auc[i] <- this_auc_higher
}
