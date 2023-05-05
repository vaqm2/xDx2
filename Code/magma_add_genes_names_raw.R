#!/usr/bin/env Rscript

library(dplyr)

args = commandArgs(trailingOnly = TRUE)

genes = read.table("/faststorage/project/xdx2/data/magma/NCBI37.3.gene.loc", 
                   header = F)
colnames(genes) = c("CODE", "CHR", "START", "END", "STRAND", "GENE")
genes = genes %>% select(CODE, GENE)
result = read.table(args[1], header = F, fill = T)
result = result %>% rename(CODE = V1)
result = inner_join(genes, result, by = c("CODE")) %>% 
    select(-CODE)

write("# VERSION = 110\n# COVAR = NSAMP MAC", args[2])

write.table(result, 
            args[2], 
            row.names = F, 
            col.names = F, 
            na = "",
            quote = F, 
            sep = " ",
            append = TRUE)