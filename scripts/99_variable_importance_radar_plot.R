# Create variable importance radar plot

setwd("/mnt/home/ryan/humpbacks/megaPrediction")

library(tidyr)
library(dplyr)
library(ggradar)
library(pals)

# Get results
fls <- list.files("./output/fitted_models_varimp/", full.names = T)
dat <- do.call(rbind, lapply(fls, read.csv, stringsAsFactors = F))

dat_sub <- dplyr::filter(dat,
                         Model == "Mr_Atlantic" |
                           Model == "Mr_EastIndian" |
                           Model == "Mr_EastPacific" |
                           Model == "Mr_Pacific" |
                           Model == "Mr_WestPacific")

plot_data <- pivot_wider(data = dat_sub,
                         names_from = "Covariate",
                         values_from = "Overall")

tolcols <- c("#0077BB",
             "#33BBEE",
             "#009988",
             "#EE7733",
             "#CC3311")

pdf("./figures/varimp_radar.pdf",
    useDingbats = F,
    width = 10, height = 6)
ggradar(plot_data,
        values.radar = c("0", "50", "100"),
        grid.min = 0,
        grid.mid = 50,
        grid.max = 100,
        gridline.mid.colour = "grey",
        group.colours = brewer.set1(5)
)        
dev.off()