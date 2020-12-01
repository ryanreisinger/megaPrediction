setwd("/mnt/home/ryan/humpbacks/megaPrediction")

## Libraries
library(ranger)
library(caret)
library(DALEX)

library(pals)

#-------------------------------------------------------------------
# Load tracking data
dat <- readRDS("./output/tracks_with_envars.RDS")

# Depending on the model, caret can only run with complete cases,
# this also means EKE, SSHA and SSHGRAD are excluded
dat <- dat[,1:19]
dat <- dat[complete.cases(dat[,9:19]), ]

# Convert response to binomial
dat[dat$real == "O", ]$real <- 1
dat[dat$real == "S", ]$real <- 0
dat$real <- as.integer(dat$real)

#-------------------------------------------------------------------
# Load data for prediction
area_grid <- readRDS("./output/envar_frame.RDS")
area_grid$EKE <- NULL
area_grid$SSHA <- NULL
area_grid$SSHGRAD <- NULL

#-------------------------------------------------------------------
# Load models
Mr_Atlantic <- readRDS("./output/fitted_models/Mr_Atlantic.RDS")
Mr_EastIndian <- readRDS("./output/fitted_models/Mr_EastIndian.RDS")
Mr_EastPacific <- readRDS("./output/fitted_models/Mr_EastPacific.RDS")
Mr_Pacific <- readRDS("./output/fitted_models/Mr_Pacific.RDS")
Mr_WestPacific <- readRDS("./output/fitted_models/Mr_WestPacific.RDS")

# Create explainers
# explainer_Atlantic <- DALEX::explain(model = Mr_Atlantic, 
#                                 data = dat[!is.na(dat$region) & dat$region == "Atlantic",9:19],
#                                 y = dat[!is.na(dat$region) & dat$region == "Atlantic", ]$real, 
#                                 label = "Mr_Atlantic")
# 
# explainer_EastIndian <- DALEX::explain(model = Mr_EastIndian, 
#                                        data = dat[!is.na(dat$region) & dat$region == "EastIndian",9:19],
#                                        y = dat[!is.na(dat$region) & dat$region == "EastIndian", ]$real, 
#                                         label = "Mr_EastIndian")
# 
# explainer_EastPacific <- DALEX::explain(model = Mr_EastPacific, 
#                                         data = dat[!is.na(dat$region) & dat$region == "EastPacific",9:19],
#                                         y = dat[!is.na(dat$region) & dat$region == "EastPacific", ]$real, 
#                                        label = "Mr_EastPacific")
# 
# explainer_Pacific <- DALEX::explain(model = Mr_Pacific, 
#                                     data = dat[!is.na(dat$region) & dat$region == "Pacific",9:19],
#                                     y = dat[!is.na(dat$region) & dat$region == "Pacific", ]$real, 
#                                        label = "Mr_Pacific")
# 
# explainer_WestPacific <- DALEX::explain(model = Mr_WestPacific, 
#                                         data = dat[!is.na(dat$region) & dat$region == "WestPacific",9:19],
#                                         y = dat[!is.na(dat$region) & dat$region == "WestPacific", ]$real, 
#                                        label = "Mr_WestPacific")


# All data?
explainer_Atlantic <- DALEX::explain(model = Mr_Atlantic,
                                     data = dat[,9:19],
                                     y = dat$real,
                                     label = "Mr_Atlantic")

explainer_EastIndian <- DALEX::explain(model = Mr_EastIndian,
                                       data = dat[,9:19],
                                       y = dat$real,
                                       label = "Mr_EastIndian")

explainer_EastPacific <- DALEX::explain(model = Mr_EastPacific,
                                        data = dat[,9:19],
                                        y = dat$real,
                                        label = "Mr_EastPacific")

explainer_Pacific <- DALEX::explain(model = Mr_Pacific,
                                    data = dat[,9:19],
                                    y = dat$real,
                                    label = "Mr_Pacific")

explainer_WestPacific <- DALEX::explain(model = Mr_WestPacific,
                                        data = dat[,9:19],
                                        y = dat$real,
                                        label = "Mr_WestPacific")
#-------------------------------------------------------------------
# Plots for n variables

envars <- data.frame("these_vars" = c("DEPTH",
                                            "SLOPE",
                                            "SLOPEDIST",
                                            "SHELFDIST",
                                            "SST",
                                            "SSTFRONT",
                                            "SSTVAR",
                                            "ICE",
                                            "ICEDIST",
                                            "ICEVAR"),
                           "title" = c("Depth (m)",
                                       "Bottom slope (deg)",
                                       "Distance to slope (km)",
                                       "Distance to shelf (km)",
                                       "Sea surface temperature (°C)",
                                       "Sea surface temperature gradient (°C/km)",
                                       "Sea surface temperature variance (°C)",
                                       "Sea ice concentration (%)",
                                       "Distance to sea ice edge (km)",
                                       "Sea ice concentration variance (%)"))

for (i in 1:length(envars$these_vars)) {
  this_var <- envars$these_vars[i]
  print(this_var)
  
  min_lim <- min(dplyr::select(dat, this_var), na.rm = T)
  max_lim <- max(dplyr::select(dat, this_var), na.rm = T)

  # Create pdps
pdp_Atlantic <- model_profile(explainer = explainer_Atlantic, variables = this_var)
pdp_EastIndian <- model_profile(explainer = explainer_EastIndian, variables = this_var)
pdp_EastPacific <- model_profile(explainer = explainer_EastPacific, variables = this_var)
pdp_Pacific <- model_profile(explainer = explainer_Pacific, variables = this_var)
pdp_WestPacific <- model_profile(explainer = explainer_WestPacific, variables = this_var)

# Plot
if (this_var == "SST") {
p <- plot(pdp_Atlantic, pdp_EastIndian, pdp_EastPacific, pdp_Pacific, pdp_WestPacific,
     title = "", subtitle = "") +
  scale_color_manual(values = brewer.set1(5), name = "Model") +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(min_lim, max_lim)) +
  labs(x = envars$title[i], y = "p(Observed track)") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black"))

} else {
  p <- plot(pdp_Atlantic, pdp_EastIndian, pdp_EastPacific, pdp_Pacific, pdp_WestPacific,
            title = "", subtitle = "") +
    scale_color_manual(values = brewer.set1(5), name = "Model") +
    # scale_y_continuous(limits = c(0, 1)) +
    # scale_x_continuous(limits = c(min_lim, max_lim), expand = c(0,0)) +
    labs(x = envars$title[i], y = "p(Observed track)") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black")) 
}

# Save
pdf(paste0("./figures/pdp/pdp_", this_var, ".pdf"), height = 3, width = 5, useDingbats = F)
print(p)
dev.off()

}
