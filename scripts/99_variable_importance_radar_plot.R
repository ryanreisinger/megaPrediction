# Create variable importance radar plot

setwd("/mnt/home/ryan/humpbacks/megaPrediction")

library(tidyr)
library(dplyr)
library(ggradar)
library(pals)

# Get results
fls <- list.files("./output/fitted_models_varimp/", full.names = T)
dat <- do.call(rbind, lapply(fls, read.csv, stringsAsFactors = F))

# Select and format data
dat_sub <- dplyr::filter(dat,
                         Model == "Mr_Atlantic" |
                           Model == "Mr_EastIndian" |
                           Model == "Mr_EastPacific" |
                           Model == "Mr_Pacific" |
                           Model == "Mr_WestPacific")

plot_data <- pivot_wider(data = dat_sub,
                         names_from = "Covariate",
                         values_from = "Overall")

# Colours
these_cols <- c("#0077BB",
             "#33BBEE",
             "#009988",
             "#EE7733",
             "#CC3311")

these_cols <- brewer.set1(5)

# Plot
pdf("./figures/varimp_radar.pdf",
    useDingbats = F,
    width = 8, height = 5)
ggradar(plot_data,
        values.radar = c("0", "50", "100"),
        grid.min = 0,
        grid.mid = 50,
        grid.max = 100,
        gridline.mid.colour = "grey",
        group.colours = these_cols,
        base.size = 9,
        group.line.width = 1,
        group.point.size = 2,
        background.circle.colour = "white",
        gridline.min.linetype = 1,
        gridline.mid.linetype = 1,
        gridline.max.linetype = 1
)        
dev.off()

# Check average varimp
dat_group <- group_by(dat_sub, Covariate) %>% 
  summarise(., avg = mean(Overall)) %>% 
  arrange(., avg)
