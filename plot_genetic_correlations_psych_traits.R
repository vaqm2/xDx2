#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    require(dplyr)
    require(ggplot2)
})

setwd("/Users/vapp0002/Documents/xDx_2/")

rG = readxl::read_xlsx("iPSYCH2015_rG_LDSC.xlsx")
rG = rG %>% 
    mutate(rG_low = round(rG - 1.96 * SE, 2), 
           rG_high = round(rG + 1.96 * SE, 2)) %>%
    mutate(COHORT = paste0(COHORT1, "_", COHORT2)) %>% 
    select(-COHORT1, -COHORT2)
rG$rG = round(rG$rG, 2)

rG$COHORT = gsub("iPSYCH2015_iPSYCH2015", "iPSYCH", rG$COHORT)
rG$COHORT = gsub("PGC_No_iPSYCH_PGC_No_iPSYCH", 
                 "External excluding iPSYCH", 
                 rG$COHORT)
rG$COHORT = gsub("PGC_PGC", "External", rG$COHORT)
rG$COHORT = gsub("iPSYCH2015_PGC_No_iPSYCH", "iPSYCH, 
                 External excluding iPSYCH", 
                 rG$COHORT)

png("Figure1B.png", res = 300, width = 10, height = 8, units = "in")

ggplot(rG, aes(x = TRAIT1, y = TRAIT2, fill = rG)) + 
    geom_tile(color = "black") + 
    theme_classic() + 
    facet_wrap(. ~ COHORT) + 
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