#!/usr/bin/env Rscript 

library(dplyr)

args = commandArgs(trailingOnly = T)
files = list.files(path = getwd(), pattern = "*.genes.out")
file.create(args[1])

write(paste("GENE", "ZSTAT", "P", "P_ADJ", "TEST", sep = "\t"), 
      args[1], 
      append = TRUE)

for (file in files) {
    print(paste("Processing", file, "....", sep = " "))
    association = read.table(file, header = T)
    test = gsub("^iPSYCH2015_EUR_", "", file)
    test = gsub("\\..*", "", test)
    significant = association %>% 
        mutate(P_ADJ = p.adjust(P, method = c("bonferroni"))) %>%
        filter(P_ADJ < 0.05) %>%
        select(GENE, ZSTAT, P, P_ADJ) %>%
        mutate(TEST = test)
    
    if(nrow(significant) > 0) {
        write.table(significant, 
                    args[1],
                    row.names = F, 
                    col.names = F, 
                    quote = F, 
                    sep = "\t",
                    append = TRUE)
    }
}