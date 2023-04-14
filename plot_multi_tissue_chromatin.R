#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)

setwd("/Users/vapp0002/Documents/xDx_2/ldsc_seg/")

multi_tissue_chromatin = read.table("Multi_tissue_chromatin_results.txt", header = T)
multi_tissue_chromatin$Trait = gsub("_CC", "", multi_tissue_chromatin$Trait)
multi_tissue_chromatin$Trait = gsub("_", " vs ", multi_tissue_chromatin$Trait)
multi_tissue_chromatin$Type = gsub(0, "Case vs Cohort", multi_tissue_chromatin$Type)
multi_tissue_chromatin$Type = gsub(1, "Case vs Other Cases", multi_tissue_chromatin$Type)
multi_tissue_chromatin$Type = gsub(2, "Case vs Case Pairwise", multi_tissue_chromatin$Type)
multi_tissue_chromatin = multi_tissue_chromatin %>% 
    mutate(Significant = ifelse(p_fdr <= 0.05, "Yes", "No"))

threshold = multi_tissue_chromatin %>% 
    arrange(p_fdr) %>% 
    filter(p_fdr <= 0.05) %>% 
    tail(1) %>% 
    select(P) %>% 
    as.numeric()

png("Cross_Disorder_LDSC_SEG_Multi_Tissue_Chromatin.png",
    res = 300, 
    units = "in",
    width = 8,
    height = 8)

ggplot(multi_tissue_chromatin %>% filter(Trait == "xDx"),
       aes(x = Name, fill = Category, size = -log10(P), y = -log10(P))) + 
    geom_point(shape = 21, color = "black") + 
    theme_bw() +
    geom_label_repel(aes(label = ifelse(Significant == "Yes", Name, "")), 
                    size = 3,
                    min.segment.length = 0, 
                    seed = 42, 
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    arrow = arrow(length = unit(0.01, "npc")),
                    color = "black",
                    nudge_x = .15,
                    nudge_y = .5,) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(face = "bold"),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(size = "none") +     
    scale_fill_manual(values = c("red", 
                                 "blue",
                                 "green",
                                 "yellow",
                                 "orange",
                                 "cyan",
                                 "pink",
                                 "gray",
                                 "violet")) +
    xlab("Tissue") +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    geom_hline(yintercept = -log10(threshold), lty = 2)

dev.off()

png("Case_vs_Others_Cohort_LDSC_SEG_Multi_Tissue_Chromatin.png",
    res = 300, 
    units = "in",
    width = 18,
    height = 10)

ggplot(multi_tissue_chromatin %>% 
           filter(Type != "Case vs Case Pairwise" & Trait != "xDx"), 
       aes(x = Name, 
           fill = Category, 
           size = -log10(P), 
           y = -log10(P))) + 
    geom_point(shape = 21) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    guides(size = "none") + 
    scale_y_continuous(breaks = seq(0, 10, 1)) + 
    geom_hline(yintercept = -log10(threshold), lty = 2) + 
    scale_fill_manual(values = c("red", 
                                 "blue", 
                                 "green", 
                                 "yellow", 
                                 "orange", 
                                 "cyan", 
                                 "pink", 
                                 "gray", 
                                 "violet")) + 
    geom_label_repel(aes(label = ifelse(Significant == "Yes", Name, ""),
                         point.size = -log10(P)), 
                     size = 3,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     arrow = arrow(length = unit(0.01, "npc")),
                     color = "black",
                     nudge_x = .15,
                     nudge_y = .5,
                     max.time = 1, 
                     max.iter = 1e5) +
    theme(axis.line = element_blank(), 
          legend.title = element_blank(), 
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    facet_grid(Trait ~ Type, 
               scales = "free",
               space = "free") +
    xlab("Tissue")

dev.off()

png("Case_vs_Case_Pairwise_LDSC_SEG_Multi_Tissue_Chromatin.png",
    res = 300, 
    units = "in",
    width = 10,
    height = 10)

ggplot(multi_tissue_chromatin %>% 
           filter(Type == "Case vs Case Pairwise"), 
       aes(x = fct_reorder(Name, Category), 
           fill = Category, 
           size = -log10(P), 
           y = -log10(P))) + 
    geom_point(shape = 21) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    guides(size = "none") + 
    geom_hline(yintercept = -log10(threshold), lty = 2) + 
    scale_fill_manual(values = c("red", 
                                 "blue", 
                                 "green", 
                                 "yellow", 
                                 "orange", 
                                 "cyan", 
                                 "pink", 
                                 "gray", 
                                 "violet")) + 
    geom_label_repel(aes(label = ifelse(Significant == "Yes", paste0(Trait, "\n", Name), ""),
                         point.size = -log10(P)), 
                     size = 3,
                     min.segment.length = 0, 
                     hjust = 1,
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     arrow = arrow(length = unit(0.01, "npc")),
                     color = "black",
                     nudge_x = .15,
                     nudge_y = .5,
                     max.time = 1, 
                     max.iter = 1e5) +
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    xlab("Tissue") + 
    scale_y_continuous(breaks = seq(0, 12, 1))

dev.off()


significant_annotations = multi_tissue_chromatin %>% 
    filter(p_fdr <= 0.05 & Category == "CNS") %>% 
    select(Name) %>% 
    unique()

significant_annotations = inner_join(multi_tissue_chromatin, 
                                     significant_annotations, by = c("Name"))

png("Significant_Tissues_LDSC_SEG_Multi_Tissue_Chromatin.png",
    res = 300, 
    units = "in",
    width = 12,
    height = 8)

ggplot(significant_annotations, aes(y = Name, x = Trait, fill = Coefficient)) + 
    geom_tile(color = "black") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(angle = 90, face = "bold", hjust = 1)) +
    geom_text(aes(label = ifelse(p_fdr <= 0.05, "*", "")), size = 4) +
    facet_grid(Category ~ Type, scales = "free", space = "free") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    xlab("") +
    ylab("Tissue")

dev.off()