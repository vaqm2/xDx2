#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(stringr)
    library(janitor)
    library(reshape2)
    library(broom)
})

setwd("/Users/vapp0002/Documents/xDx_2/Magma/")

magma = fread("/Users/vapp0002/Documents/xDx_2/Magma/C_Magma.txt", header = T)
enriched_genes = magma %>% 
    filter(P_ADJ < 0.05) %>% 
    group_by(TEST) %>%
    mutate(N_GENES = n()) %>%
    ungroup() %>% 
    filter(N_GENES >= 20) %>% 
    select(-N_GENES)

metadata = fread("/Users/vapp0002/Documents/xDx_2/GSE25219/metadata_GSE25219.tsv",
                 header = T) %>% 
    clean_names() %>% 
    select(refinebio_accession_code, 
           characteristics_ch1_stage, 
           characteristics_ch1_age) %>%
    rename(SAMPLE = refinebio_accession_code,
           STAGE = characteristics_ch1_stage,
           AGE = characteristics_ch1_age)

# metadata$STAGE = factor(metadata$STAGE)

gse25219 = fread("/Users/vapp0002/Documents/xDx_2/GSE25219/GSE25219.tsv", 
                 header = T,
                 sep = "\t")

gse25219 = gse25219 %>% 
    melt(id.vars = c("Gene")) %>% 
    rename(SAMPLE = variable) %>%
    rename(GENE = Gene) %>%
    rename(normalized_expression = value) # %>% 
#    group_by(SAMPLE) %>%
#    mutate(normalized_expression_scaled = scale(normalized_expression, 
#                                                center = T, 
#                                                scale = F)) %>%
#    ungroup()

gse25219_genes_of_interest = inner_join(gse25219, 
                                        enriched_genes, 
                                        by = c("GENE"), 
                                        relationship = "many-to-many") %>% 
    group_by(TEST, SAMPLE) %>% 
    summarise(MeanExpression = mean(normalized_expression), 
              .groups = "drop")

gse25219_genes_of_interest = inner_join(gse25219_genes_of_interest, metadata, 
                                        by = c("SAMPLE"),
                                        relationship = "many-to-many")

png("Developmental_Trajectories.png", 
    width = 10, 
    height = 8, 
    units = "in", 
    res = 300)

ggplot(gse25219_genes_of_interest, aes(x = STAGE, y = MeanExpression, color = TEST)) + 
    geom_smooth(method = "loess") +
    theme_bw() +
    scale_x_continuous(breaks = seq(1, 20, by = 1)) +
    scale_y_continuous(breaks = seq(0, 10, by = c(0.1))) +
    scale_color_manual(values = c("red", "blue", "green", "black", "yellow", 
                                  "orange", "cyan", "brown", "violet")) + 
    ylab("Normalized Expression") +
    theme(legend.title = element_blank()) + 
    geom_vline(xintercept = 8, lty = 2, color = "red") +
    geom_hline(yintercept = 0)

dev.off()