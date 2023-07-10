#!/usr/bin/env Rscript

library(dplyr)
library(CMplot)
library(vroom)
library(purrr)

args = commandArgs(trailingOnly = TRUE)
files = vroom(args[1])
list_of_assoc = map(files, function(file = path, new_col = trait) {
    this_assoc = vrooom(path) %>% 
        select(SNP, CHR, BP, P) %>%
        mutate(Trait = trait)
    return(this_assoc)
})

assoc = bind_rows(list_of_assoc)

CMplot(assoc, 
       type = "p",
       plot.type = c("c", "q"),
       col = c("black", "red", "blue", "green", "brown", "cyan", "orange"),
       outward = TRUE,
       threshold = 5e-8,
       file = "tiff", 
       file.name = paste0(args[2], ".tiff"),
       dpi = 300,
       file.output = TRUE,
       verbose = TRUE,
       width = 15,
       height = 15)