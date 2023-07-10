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

assoc_0.05 = bind_rows(list_of_assoc) %>%
    filter(P < 0.05) %>%
    pivot_wider(names_from = Trait, values_from = P)

assoc = bind_rows(list_of_assoc) %>%
    pivot_wider(names_from = Trait, values_from = P)

png(paste0(args[2], "_0.05_Manhattan.png"), 
    res = 300,
    width = 15,
    height = 15,
    units = "in")

CMplot(assoc, 
       type = "p",
       plot.type = "c",
       col = c("grey60", "grey30"),
       outward = FALSE,
       r = 2,
       cir.chr.h = 2,
       multracks = TRUE,
       threshold = 5e-8,
       threshold.col = "blue",
       trait.legend.ncol = 7,
       signal.line = 1,
       pch = 21,
       signal.pch = 24,
       amplify = FALSE,
       signal.col = "red",
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       dpi = 300,
       file.output = TRUE,
       conf.int = TRUE,
       width = 15,
       height = 15)

dev.off()

png(paste0(args[2], "_0.05_Manhattan.png"), 
    res = 300,
    width = 15,
    height = 15,
    units = "in")

CMplot(assoc, 
       type = "p",
       plot.type = "m",
       col = c("grey60", "grey30"),
       multracks = TRUE,
       threshold = c(1e-6, 5e-8),
       threshold.col = c("green", "blue"),
       trait.legend.ncol = 7,
       signal.line = 1,
       pch = 21,
       signal.pch = 24,
       amplify = FALSE,
       signal.col = "red",
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       dpi = 300,
       file.output = TRUE,
       conf.int = TRUE,
       width = 15,
       height = 15)

dev.off()

png(paste0(args[2], "_0.05_QQ.png"), 
    res = 300,
    width = 10,
    height = 10,
    units = "in")

CMplot(assoc, 
       type = "p",
       plot.type = "q",
       signal.pch = c(21:27),
       signal.cex = 1.2,
       signal.col = "black",
       col = c("red", "blue", "green", "yellow", "cyan", "orange", "grey"),
       threshold = 5e-8,
       trait.legend.ncol = 7,
       pch = 21,
       trait.legend.pos = "right",
       ylab=expression(-log[10](italic(p))),
       multracks = TRUE,
       dpi = 300,
       file.output = TRUE,
       conf.int = TRUE,
       width = 8,
       height = 8)

dev.off()