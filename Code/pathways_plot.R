#!/usr/bin/env Rscript

library(janitor)
library(ggplot2)
library(ggpattern)

setwd("/Users/vapp0002/Documents/xDx_2/Magma/Pathways/")

pathways = read.table("GO_sets.txt", 
                      header = T)
pathways$SOURCE = sub("^(.*?)_.*", "\\1", pathways$SET)
pathways$SOURCE = gsub("^GO", "", pathways$SOURCE)
pathways$SET = sub("^[^_]*_", "", pathways$SET)
pathways$SET = make_clean_names(pathways$SET, allow_dupes = TRUE)
pathways$SET = gsub("_", " ", pathways$SET)
pathways = pathways %>% mutate(TYPE = ifelse(str_detect(TEST, "_CC"),
                                             "Case vs. Other Cases",
                                             "Case vs. Cohort"))

pathways_top10 = pathways %>% 
    group_by(TEST) %>% 
    arrange(P) %>% 
    slice(1:10) %>%
    ungroup() %>%
    select(SET, TEST) %>%
    unique()

top_pathways = inner_join(pathways, pathways_top10, by = c("SET", "TEST"))

top_pathways$TEST = gsub("_CC", "", top_pathways$TEST)

png("GO_Pathways_C_Magma.png",
    res = 300,
    width = 18,
    height = 15,
    units = "in")

ggplot(top_pathways, aes(x = -log10(P), y = SET, fill = SOURCE)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(TEST ~ TYPE, scales = "free") +
    geom_vline(xintercept = -log10(0.05/length(unique(pathways$SET))), lty = 2) +
    theme_bw() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_fill_manual(values = c("red", "blue", "green")) +
    ylab("") +
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom")

dev.off()