#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(readxl)
library(janitor)
library(clusterProfiler)
library(org.Hs.eg.db)

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

iso_twas_adult = read_all_sheets("iPSYCH_CC_Finemapped_isoTWAS_Adult.xlsx")
iso_twas_dev   = read_all_sheets("iPSYCH_CC_Finemapped_isoTWAS_Developmental.xlsx")

gene_list = iso_twas_adult %>% 
    filter(trait == "ADHD_CC") %>% 
    arrange(desc(z)) %>%
    pull(z)
names(gene_list) = iso_twas_adult %>% 
    filter(trait == "ADHD_CC") %>% 
    arrange(desc(z)) %>%
    pull(gene)

gse = gseGO(geneList = gene_list,
            keyType = "ENSEMBL",
#            nperm = 10000,
            OrgDb = org.Hs.eg.db,
            ont = "ALL",
            pAdjustMethod = "fdr",
            pvalueCutoff = 0.05,
            minGSSize = 10,
            maxGSSize = 500)

as.data.frame(gse)