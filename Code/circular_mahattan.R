#!/usr/bin/env Rscript

library(tidyverse)
library(CMplot)

args = commandArgs(trailingOnly = TRUE)
files = vroom(args[1])
list_of_assoc = map2(files %>% pull(Path), 
                     files %>% pull(Trait),
                     function(path, trait) {
    this_assoc = vroom(path) %>% 
        select(SNP, CHR, BP, P) %>%
        mutate(Trait = trait)
    return(this_assoc)
})

assoc = bind_rows(list_of_assoc) %>%
    pivot_wider(names_from = Trait, values_from = P)

CMplot(assoc, 
       type = "p",
       plot.type = c("c", "q"),
       col = c("black", "red", "blue", "green", "brown", "cyan", "orange"),
       outward = TRUE,
       threshold = 5e-8,
       file = "tiff", 
       file.name = args[2],
       dpi = 300,
       file.output = TRUE,
       verbose = TRUE,
       width = 15,
       height = 15)