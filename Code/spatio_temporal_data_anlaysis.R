#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(stringr)
    library(janitor)
    library(reshape2)
})

metadata = fread("/Users/vapp0002/Downloads/GSE25219/metadata_GSE25219.tsv",
                 header = T) %>% 
    clean_names() %>% select(refinebio_accession_code, 
                             contains("characteristics"),
                             -characteristics_ch1_ph,
                             -characteristics_ch1_postmortem_interval,
                             -characteristics_ch1_rna_integrity_number)

gse25219 = fread("/Users/vapp0002/Downloads/GSE25219/GSE25219.tsv", 
                 header = T,
                 sep = "\t")

magma = fread("/Users/vapp0002/Documents/xDx_2/Magma/C_Magma.txt", header = T)

ggplot(metadata, aes(x = as.factor(characteristics_ch1_stage))) + 
    geom_bar() + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("")

age_stage_df = metadata %>% 
    select(characteristics_ch1_age, characteristics_ch1_stage) %>% 
    table() %>%
    as.data.frame() %>% 
    filter(Freq > 0)

age_order = metadata %>% 
    select(characteristics_ch1_age, characteristics_ch1_stage) %>%
    unique() %>% 
    arrange(characteristics_ch1_stage, characteristics_ch1_age) %>% 
    select(characteristics_ch1_age)

ggplot(age_stage_df, aes(y = factor(characteristics_ch1_stage), 
                         x = characteristics_ch1_age, 
                         fill = Freq)) + 
    geom_tile(color = "black") + 
    theme_bw() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    geom_text(aes(label = ifelse(Freq > 0, Freq, ""))) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_x_discrete(limits = age_order$characteristics_ch1_age) +
    xlab("Age") + 
    ylab("Stage")

gse25219 = gse25219 %>% 
    melt(id.vars = c("Gene")) %>% 
    rename(refinebio_accession_code = variable) %>%
    rename(normalized_expression = value)

gse25219 = inner_join(gse25219, metadata, by = c("refinebio_accession_code"))

head(gse25219)