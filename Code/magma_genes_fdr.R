#!/usr/bin/env Rscript 

library(dplyr)

setwd("/Users/vapp0002/Documents/xDx_2/Magma/C_Magma_Genes/")
args = commandArgs(trailingOnly = T)
files = list.files(path = getwd(), pattern = "*hgnc.out")

genes = data.frame(GENE = character(),
                   Z = numeric(),
                   P = numeric(),
                   TEST = character())

for (file in files) {
    print(paste("Processing", file, "....", sep = " "))
    association = read.table(file, header = T)
    test = gsub("^iPSYCH2015_EUR_", "", file)
    test = gsub("\\..*", "", test)
    association = association %>% 
        select(GENE, ZSTAT, P) %>% 
        rename(Z = ZSTAT) %>% 
        mutate(TEST = test)
    genes = rbind(genes, association)
}

genes = genes %>% 
    mutate(P_ADJ = p.adjust(P, method = c("fdr")))
genes_significant = genes %>%    
    filter(P_ADJ <= 0.05) %>% 
    select(GENE) %>% 
    unique()

if(nrow(genes_significant) > 0) {
    genes_significant = inner_join(genes, genes_significant, 
                                   by = c("GENE"),
                                   relationship = "many-to-many")
    write.table(genes_significant, 
                "GenesFDR.txt", 
                row.names = F, 
                sep = "\t", 
                quote = F)
} else {
    print("No significant genes after multiple-testing correction!")
}