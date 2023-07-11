#!/usr/bin/env Rscript

library(vroom)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

significant_genes = function(x_tbl, fdr_threshold = 0.05) {
    genes_pass_fdr = x_tbl %>% 
        mutate(P_FDR = p.adjust(P, method = "fdr")) %>%
        filter(P_FDR <= fdr_threshold) %>% 
        pull(GENE)
    return(genes_pass_fdr)
}

args = commandArgs(trailingOnly = TRUE)
test_set = vroom(args[1], show_col_types = FALSE) %>% significant_genes()
# background = vroom(args[2], show_col_types =FALSE) %>% significant_genes()

setwd("/Users/vapp0002/Documents/xDx_2/Magma/C_Magma_Genes/")

go = enrichGO(gene = test_set,
              ont = "ALL",
              OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL",
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 0.01,
              qvalueCutoff = 0.05,  
              pAdjustMethod = "fdr")
              # universe = background)

as.data.frame(go)