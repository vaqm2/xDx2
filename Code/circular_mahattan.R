#!/usr/bin/env Rscript

library(dplyr)
library(CMplot)
library(vroom)
library(tidyr)
library(purrr)

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
       col = c("grey30", "grey60"),
       outward = FALSE,
       threshold = c(1e-6, 5e-8),
       threshold.col = c("blue", "red"),
       trait.legend.ncol = 7,
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       file = "tiff", 
       file.name = args[2],
       dpi = 300,
       file.output = TRUE,
       verbose = TRUE,
       conf.int = TRUE,
       width = 15,
       height = 15)