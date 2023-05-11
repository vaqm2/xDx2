#!/usr/bin/env Rscript 

suppressPackageStartupMessages({
library(dplyr)
library(ggplot2)})

setwd("/Users/vapp0002/Documents/xDx_2/")

rG = readxl::read_xlsx("xDx_rG_comparison.xlsx", col_names = TRUE)
rG = rG %>% filter(METHOD == "S-PCGC" & 
                       COVAR_REGRESS == "PCs" &
                       PHENO1 != "xDx" & PHENO2 != "xDx") %>%
    mutate(rG_high = round(rG + 1.96 * SE, 2), rG_low = round(rG - 1.96 * SE), 2)

png("Publication/Figure1A.png", res = 300, width = 5, height = 4, units = "in")

ggplot(rG, aes(x = PHENO1, y = PHENO2, fill = rG)) + 
    geom_tile(color = "black") + 
    theme_classic() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") + 
    xlab("") + 
    ylab("") + 
    geom_text(aes(label = paste0(rG, "\n", rG_low, " - ", rG_high)), 
              size = 2.75) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"),
          legend.text = element_text(face = "bold")) +
    guides(fill = guide_colorbar(frame.colour = "black"))

dev.off()