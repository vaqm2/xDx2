#!/usr/bin/env Rscript 

library(dplyr)

args = commandArgs(trailingOnly = T)
files = list.files(path = getwd(), pattern = "*.gsa.out")

gene_sets = data.frame(SET = character(),
                       Z = numeric(),
                       P = numeric(),
                       TEST = character())

for (file in files) {
    print(paste("Processing", file, "....", sep = " "))
    association = read.table(file, skip = 3, header = T)
    test = gsub("^iPSYCH2015_EUR_", "", file)
    test = gsub("\\..*", "", test)
    association = association %>% 
        filter(NGENES >= 10 & NGENES <= 500) %>%
        mutate(TEST = test) %>% 
        mutate(Z = BETA_STD/SE) %>%
        select(-TYPE, -BETA, -NGENES, -BETA_STD, -SE)
    if(FULL_NAME %in% colnames(association)) {
        association = association %>% 
            rename(SET = FULL_NAME) %>% 
            select(-VARIABLE)
    } else {
        association = association %>% 
            rename(SET = VARIABLE)
    }
    gene_sets = rbind(gene_sets, association)
}

gene_sets %>% mutate(P_ADJ = p.adjust(P, method = c("fdr")))

write.table(gene_sets, "GeneSetsFDR.txt", row.names = F, sep = "\t", quote = F)