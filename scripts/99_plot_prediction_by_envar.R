library(dplyr)
library(ggplot2)
library(tidyr)

setwd("C:/Users/Ryan Reisinger/Documents/Academic/UCSC/Work/Analysis/megaPrediction")

# Option 1: plotting predictions against covariates
# This is NOT a partial dependence plot
dat <- readRDS("./output/predictions/M5_grid.RDS")

dat <- dplyr::select(dat,
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
                     Mr_Pacific_p,
                     Mr_EastPacific_p,
                     Mr_WestPacific_p,
                     Mr_EastPacific_p,
                     Mr_Atlantic_p)

dat_long <- tidyr::pivot_longer(data = dat,
                                cols = DEPTH:ICEVAR,
                                names_to = "Covariate",
                                values_to = "Value")

dat_long <- tidyr::pivot_longer(data = dat_long,
                    cols = Mr_Pacific_p:Mr_Atlantic_p,
                    names_to = "Model",
                    values_to = "p")

tiff("./figures/response_curves.tiff",
     height = 8, width = 12, res = 300, units = "in")
ggplot(data = dat_long, aes(y = p, x = Value)) +
  geom_bin2d() +
  geom_smooth() +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  facet_grid(rows = vars(Model), cols = vars(Covariate), scales = "free")
dev.off()