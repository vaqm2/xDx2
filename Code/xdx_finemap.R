#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(readxl)
library(janitor)
library(cowplot)
library(UpSetR)

setwd("/Users/vapp0002/Documents/xDx_2/finemap/")

read_all_sheets = function(xl_file, bind = TRUE) {
    sheets = excel_sheets(xl_file)
    sheet_tbl = lapply(sheets, 
                       function(sheet_name) { 
                           read_xlsx(xl_file, sheet = sheet_name) %>%
                               clean_names()
                           }
                       )
    names(sheet_tbl) = sheets
    if(bind) sheet_tbl = bind_rows(sheet_tbl) 
    return(sheet_tbl)
}

label_types = function(x_df) {
    x_df %>% mutate(Type = str_count(Trait, "_")) %>%
        mutate(Type = gsub("0", " vs. Cohort", Type)) %>%
        mutate(Type = gsub("1", " vs. Other Cases", Type)) %>%
        mutate(Type = gsub("2", "", Type)) %>%
        mutate(Trait = gsub("_CC", "", Trait)) %>%
        mutate(Trait = gsub("_", " vs. ", Trait))
}

coloc          = read_all_sheets("ipsych_xGWAS_coloc.xlsx", bind = FALSE)
iso_twas_adult = read_all_sheets("iPSYCH_CC_Finemapped_isoTWAS_Adult.xlsx")
iso_twas_dev   = read_all_sheets("iPSYCH_CC_Finemapped_isoTWAS_Developmental.xlsx")
twas_adult     = read_all_sheets("iPSYCH_CC_Finemapped_TWAS_Adult.xlsx")
twas_dev       = read_all_sheets("iPSYCH_CC_Finemapped_TWAS_Developmental.xlsx")

coloc = coloc %>% 
    lapply(function(df) df %>% select(snp_id, clpp, gwas, contains("gene_name"))) %>%
    bind_rows(.id = "id") %>% 
    mutate(GENE = ifelse(!is.na(gene_name), gene_name, leafviz_annot_gene_name)) %>% 
    select(-gene_name, -leafviz_annot_gene_name)
iso_twas = rbind(iso_twas_adult, iso_twas_dev)
twas = rbind(twas_adult, twas_dev)

coloc_genes = coloc %>% 
    select(GENE, gwas, id, clpp) %>% 
    rename(Gene = GENE, Trait = gwas, Method = id) %>%
    label_types() %>%
    na.omit() %>%
    filter(Gene != "NA")
twas_genes = twas %>% 
    select(hgnc, z, trait, tissue) %>% 
    rename(Gene = hgnc, Trait = trait, Method = tissue) %>%
    label_types() %>%
    na.omit() %>%
    filter(Gene != "NA")
iso_twas_genes = iso_twas %>% 
    select(hgnc, z, trait, tissue) %>% 
    rename(Gene = hgnc, Trait = trait, Method = tissue) %>%
    label_types() %>%
    na.omit() %>%
    filter(Gene != "NA")

# Plotting

png("xDx_colocalization.png", res = 300, width = 12, height = 12, units = "in")

coloc_genes %>% 
    ggplot(aes(y = Gene, x = Method, fill = -log2(clpp))) +
    geom_tile(color = "black") + 
    theme_bw() +
    facet_grid(. ~ paste0(Trait, Type), scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          strip.text = element_text(angle = 90, hjust = 1),
          strip.background = element_rect(fill = "white")) +
    xlab("") + ylab("") +
    scale_fill_gradient(low = "blue", high = "red")

dev.off()

png("xDx_TWAS.png", res = 300, width = 12, height = 12, units = "in")

twas_genes %>% 
    ggplot(aes(y = Gene, x = Method, fill = z)) +
    geom_tile(color = "black") + 
    theme_bw() +
    facet_grid(. ~ paste0(Trait, Type), scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          strip.text = element_text(angle = 90, hjust = 1),
          strip.background = element_rect(fill = "white")) +
    xlab("") + ylab("") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")

dev.off()

iso_twas_genes %>% 
    ggplot(aes(x = Gene, y = paste0(Trait, Type), fill = z)) +
    geom_tile(color = "black") + 
    theme_bw() +
    facet_grid(Method ~ ., scales = "free_y", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          strip.background = element_rect(fill = "white")) +
    xlab("") + ylab("") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")