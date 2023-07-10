#!/usr/bin/env Rscript

library(dplyr)
library(CMplot)
library(vroom)
library(tidyr)
library(purrr)

args = commandArgs(trailingOnly = TRUE)
files = vroom(args[1], show_col_types = FALSE)
list_of_assoc = map2(files %>% pull(Path), 
                     files %>% pull(Trait),
                     function(path, trait) {
    this_assoc = vroom(path, show_col_types = FALSE) %>% 
        select(SNP, CHR, BP, P) %>%
        mutate(Trait = trait)
    return(this_assoc)
})

assoc_0.05 = bind_rows(list_of_assoc) %>%
    filter(P < 0.05) %>%
    pivot_wider(names_from = Trait, values_from = P)

assoc = bind_rows(list_of_assoc) %>%
    pivot_wider(names_from = Trait, values_from = P)

CMplot(assoc_0.05, 
       H = 2,
       type = "p",
       plot.type = "c",
       col = c("grey60", "grey30"),
       outward = FALSE,
       r = 2,
       cir.chr.h = 2,
       threshold = 5e-8,
       threshold.col = "blue",
       trait.legend.ncol = 7,
       signal.line = 1,
       pch = 21,
       amplify = TRUE,
       signal.cex = 1.25,
       signal.col = "red",
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       dpi = 300,
       file.output = TRUE,
       conf.int = TRUE,
       width = 15,
       height = 15,
       file = "tiff",
       cir.legend = TRUE,
       chr.border = TRUE,
       file.name = paste0(args[2], "_0.05_Manhattan"))

CMplot(assoc, 
       H = 2,
       type = "p",
       plot.type = "m",
       col = c("grey60", "grey30"),
       threshold = c(1e-6, 5e-8),
       threshold.col = c("green", "blue"),
       trait.legend.ncol = 7,
       signal.line = 1,
       pch = 21,
       signal.cex = 1.25,
       signal.col = "red",
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       dpi = 300,
       file.output = TRUE,
       amplify = TRUE,
       conf.int = TRUE,
       width = 15,
       height = 15,
       file = "tiff",
       cir.legend = TRUE,
       chr.border = TRUE,
       file.name = paste0(args[2], "_Manhattan"))

CMplot(assoc, 
       type = "p",
       plot.type = "q",
       signal.pch = c(21:27),
       signal.cex = 1.2,
       signal.col = "black",
       col = c("red", "blue", "green", "brown", "cyan", "orange", "grey"),
       threshold = 5e-8,
       trait.legend.ncol = 7,
       pch = 21,
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       dpi = 300,
       file.output = TRUE,
       conf.int = TRUE,
       width = 8,
       height = 8,
       file = "tiff",
       file.name = paste0(args[2], "_qq"))