#!/usr/bin/env Rscript

args = commandArgs(trailingOnly =  TRUE)

library(dplyr)

load("/faststorage/project/xdx2/data/magma/H-Magma/geneAnno_allgenes.rda")
geneAnno1 = geneAnno1 %>% 
    filter(gene_biotype == "protein_coding" & hgnc_symbol != "") %>% 
    select(ensembl_gene_id, hgnc_symbol)

genes_to_test = read.table(args[1], header = F, fill = T)
genes_to_test = genes_to_test %>% rename(ensembl_gene_id = V1)

genes_to_test = inner_join(geneAnno1, genes_to_test) %>% 
    select(-ensembl_gene_id)

write.table(genes_to_test, 
            args[2], 
            row.names = F, 
            col.names = F,
            quote = F, 
            sep = "\t")