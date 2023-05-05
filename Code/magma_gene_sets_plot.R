#!/usr/bin/env Rscript 

library(dplyr)
library(ggplot2)
library(scales)

args = commandArgs(trailingOnly = TRUE)

prefix = "iPSYCH2015_EUR_"
suffix = "_C5_GO.gsa.out"

xdx = read.table(paste0(prefix, "xDx", suffix), header = T)
xdx = xdx %>%
    arrange(P) %>%
    head(20) %>%
    select(VARIABLE, FULL_NAME, BETA, BETA_STD, P) %>%
    rename(GO = VARIABLE) %>%
    mutate(TRAIT = "xDx") %>% 
    mutate(GWAS = "Case vs Cohort") %>%
    mutate(BETA_HIGH = BETA + 1.96 * BETA_STD) %>%
    mutate(BETA_LOW = BETA - 1.96 * BETA_STD)
xdx$GO = gsub("_.*", "", xdx$GO)
xdx$GO = gsub("^GO", "", xdx$GO)
xdx$FULL_NAME = substr(xdx$FULL_NAME, 6, nchar(xdx$FULL_NAME))
xdx$FULL_NAME = gsub("_", " ", xdx$FULL_NAME)

case_case = xdx[FALSE,]
pairwise  = xdx[FALSE,]

for (trait in c("ADHD", "ANO", "AUT", "BIP", "MDD", "SCZ")) {
    file      = read.table(paste0(prefix, trait, suffix), header = T)
    file_cc   = read.table(paste0(prefix, trait, "_CC", suffix), header = T)
    top_5    = file %>% arrange(P) %>% head(5) %>% select(FULL_NAME)
    top_5_cc = file_cc %>% arrange(P) %>% head(5) %>% select(FULL_NAME)
    to_subset = rbind(top_5, top_5_cc) %>% unique()
    file      = inner_join(to_subset, file, by = c("FULL_NAME")) %>%
        select(VARIABLE, FULL_NAME, BETA, BETA_STD, P) %>%
        rename(GO = VARIABLE) %>%
        mutate(GWAS = "Case vs Cohort") %>%
        mutate(BETA_HIGH = BETA + 1.96 * BETA_STD) %>%
        mutate(BETA_LOW = BETA - 1.96 * BETA_STD)
    file_cc   = inner_join(to_subset, file_cc, by = c("FULL_NAME")) %>%
        select(VARIABLE, FULL_NAME, BETA, BETA_STD, P) %>%
        rename(GO = VARIABLE) %>%
        mutate(GWAS = "Case vs Other Cases") %>%
        mutate(BETA_HIGH = BETA + 1.96 * BETA_STD) %>%
        mutate(BETA_LOW = BETA - 1.96 * BETA_STD)
    merged    = rbind(file, file_cc) %>% mutate(TRAIT = trait)
    case_case = rbind(case_case, merged)
}

case_case$GO = gsub("_.*", "", case_case$GO)
case_case$GO = gsub("^GO", "", case_case$GO)
case_case$FULL_NAME = substr(case_case$FULL_NAME, 6, nchar(case_case$FULL_NAME))
case_case$FULL_NAME = gsub("_", " ", case_case$FULL_NAME)


for (trait in c("ADHD_AUT", "ADHD_ANO", "ADHD_BIP", "ADHD_MDD", "ADHD_SCZ",
                "ANO_AUT", "ANO_BIP", "ANO_MDD", "ANO_SCZ",
                "AUT_BIP", "AUT_MDD", "AUT_SCZ",
                "BIP_MDD", "BIP_SCZ",
                "MDD_SCZ")) {
    file = read.table(paste0(prefix, trait, "_CC", suffix), header = T)
    file = file %>%
        arrange(P) %>%
        head(5) %>%
        select(VARIABLE, FULL_NAME, BETA, BETA_STD, P) %>%
        rename(GO = VARIABLE) %>%
        mutate(TRAIT = trait) %>% 
        mutate(GWAS = "Case vs Case Pairwise") %>%
        mutate(BETA_HIGH = BETA + 1.96 * BETA_STD) %>%
        mutate(BETA_LOW = BETA - 1.96 * BETA_STD)
    pairwise = rbind(pairwise, file)
}

pairwise$GO = gsub("_.*", "", pairwise$GO)
pairwise$GO = gsub("^GO", "", pairwise$GO)
pairwise$FULL_NAME = substr(pairwise$FULL_NAME, 6, nchar(pairwise$FULL_NAME))
pairwise$FULL_NAME = gsub("_", " ", pairwise$FULL_NAME)

png(paste0(args[1], "_xDx.png"), res = 300, width = 8, height = 8, units = "in")

ggplot(xdx, aes(y = FULL_NAME, x = BETA, shape = GO, fill = -log10(P))) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = BETA - 1.96 * BETA_STD, 
                       xmax = BETA + 1.96 * BETA_STD),
                   height = 0.01, linewidth = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    ylab("") +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_gradient(low = "blue", high = "red")

dev.off()

png(paste0(args[1], "_Case_Case.png"), res = 300, width = 15, height = 15, units = "in")

ggplot(case_case, aes(y = FULL_NAME, 
                      x = BETA, 
                      shape = GO,
                      fill = -log10(P))) + 
    geom_point(size = 3, position = position_dodge(width = 0.9)) +
    geom_errorbarh(aes(xmin = BETA - 1.96 * BETA_STD, 
                       xmax = BETA + 1.96 * BETA_STD),
                   height = 0.01,
                   position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(size = 8, 
                                     hjust = 1, 
                                     face = "bold"),
          legend.position = "bottom") +
    scale_y_discrete(labels = label_wrap(40), guide = guide_axis(n.dodge = 2)) +
    facet_grid(TRAIT ~ GWAS, scales = "free", space = "free") + 
    scale_fill_gradient(high = "red", low = "blue") +
    ylab("") +
    scale_shape_manual(values = c(21, 22, 23)) +
    geom_vline(xintercept = 0, lty = 2)

dev.off()

png(paste0(args[1], "_Pairwise.png"), res = 300, width = 15, height = 15, units = "in")

ggplot(pairwise, aes(y = FULL_NAME, x = BETA, shape = GO, fill = -log10(P))) + 
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = BETA - 1.96 * BETA_STD, 
                       xmax = BETA + 1.96 * BETA_STD),
                   height = 0.01) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(size = 8, 
                                     hjust = 1, 
                                     face = "bold"),
          legend.position = "bottom") +
    scale_y_discrete(labels = label_wrap(30)) +
    ylab("") + 
    facet_wrap(TRAIT ~ ., scales = "free") +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_gradient(low = "blue", high = "red")

dev.off()