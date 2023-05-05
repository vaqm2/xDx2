#!/usr/bin/env Rscript

library(dplyr)

args = commandArgs(trailingOnly = TRUE)

genes = read.table("/faststorage/project/xdx2/data/magma/NCBI37.3.gene.loc", 
                   header = F)
colnames(genes) = c("CODE", "CHR", "START", "END", "STRAND", "GENE")
genes = genes %>% select(CODE, GENE)
result = read.table(args[1], header = T)
result = result %>% rename(CODE = GENE)
result = inner_join(genes, result, by = c("CODE")) %>% 
    select(-CODE) %>% 
    arrange(desc(abs(ZSTAT)))

write.table(result, 
            args[2], 
            row.names = F, 
            col.names = T, 
            quote = F, 
            sep = "\t")