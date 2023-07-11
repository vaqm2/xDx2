#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(cowplot)))
library(readxl)
library(ggpubr)
library(tidyverse)

setwd("/Users/vapp0002/Documents/xDx_2/PsychEncode/")
magma_path = "/Users/vapp0002/Documents/xDx_2/Magma/C_Magma_Genes/"

datMeta_subject <- openxlsx::read.xlsx("http://development.psychencode.org/files/raw_data/mRNA-seq_Sample%20metadata.xlsx",
                                       startRow = 2,)


datMeta_sample <- openxlsx::read.xlsx("http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_QC.xlsx",
                                      startRow = 3)
datMeta_sample <- slice(datMeta_sample, -1) %>% 
    select(2:3) %>% 
    unite(ID, c(Braincode, Regioncode), sep = ".", remove = FALSE)

datMeta <- left_join(datMeta_sample, datMeta_subject, by = "Braincode")  # 607 samples

datMeta <- datMeta %>% 
    mutate(Region = ifelse(Regioncode %in% c("DFC", "MFC", "PC", "VFC", "M1C", "OFC", "FC"), "Frontal Cortex",
                           ifelse(Regioncode %in% c("A1C", "ITC", "STC", "TC"), "Temporal Cortex",
                                  ifelse(Regioncode %in% c("IPC", "S1C", "PC"), "Parietal Cortex",
                                         ifelse(Regioncode %in% c("OC", "V1C"), "Visual Cortex",
                                                ifelse(Regioncode %in% c("CB", "CBC", "URL"), "Cerebellum",
                                                       ifelse(Regioncode %in% c("DTH", "MD", "DIE"), "Thalamus",
                                                              ifelse(Regioncode %in% c("CGE", "LGE" ,"MGE", "STR"), "Striatum",
                                                                     ifelse(Regioncode == "HIP", "Hippocampus", 
                                                                            ifelse(Regioncode == "AMY","Amygdala", NA))))))))))

datMeta$Region2 = datMeta$Region
datMeta$Region2[grep("Cortex", datMeta$Region2)] = "Cortex"
datMeta <- mutate(datMeta, Period = ifelse(Days < 40 * 7, "Prenatal", "Postnatal"))
datMeta$Period = factor(datMeta$Period, levels = c("Prenatal", "Postnatal"))
datMeta %>% head()

expr = as.data.frame(fread("mRNA_Expression.txt"))

expr = expr %>% pivot_longer(names_to = "ID", values_to = "RPKM", cols = !"Geneid")

expr_data = inner_join(datMeta, expr, by = "ID")
expr_data$Geneid <- gsub("^[^|]*\\|", "", expr_data$Geneid)
head(expr_data)

genes = data.frame(GENE = character(),
                   Z = numeric(),
                   P = numeric(),
                   TEST = character())
files = list.files(path = "/Users/vapp0002/Documents/xDx_2/Magma/C_Magma_Genes/", 
                   pattern = "*.hgnc.out")

for (file in files) {
    print(paste("Processing", file, "....", sep = " "))
    association = read.table(paste0(magma_path, file), header = T)
    test = gsub("^iPSYCH2015_EUR_", "", file)
    test = gsub("\\..*", "", test)
    association = association %>% 
        mutate(TEST = test) %>% 
        mutate(P_FDR = p.adjust(P, method = "fdr")) %>%
        filter(P_FDR <= 0.05) %>%
        select(GENE, TEST)
    if(nrow(association) >= 100) {
        genes = rbind(genes, association)
    }
}

enriched_genes = genes %>% 
    mutate(TYPE = str_count(TEST, "_")) %>%
    filter(TYPE < 2)

expr_data_magma = inner_join(expr_data, enriched_genes, 
                             by = join_by(Geneid == GENE),
                             relationship = "many-to-many")

mean_expr_by_test = expr_data_magma %>% 
    group_by(TEST, ID) %>% 
    summarise(Mean_RPKM = mean(RPKM), 
              .groups = "drop") %>% 
    mutate(TYPE = ifelse(str_detect(TEST, "_"), 
                         "Case vs. Other Cases", 
                         "Case vs. Cohort"))

png("Pre_vs_Postnatal_Expression_min_100.png",
    res = 300,
    width = 15,
    height = 10,
    units = "in")

ggplot(na.omit(expr_data_magma), aes(x = TEST, y = log(0.1 + RPKM), fill = Period)) + 
    geom_boxplot(position = position_dodge(width = 0.9)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank()) + 
    facet_wrap(Region ~ .) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    xlab("") +
    scale_fill_manual(values = c("red", "lightblue"))

dev.off()
    
mean_expr_by_test = inner_join(mean_expr_by_test, datMeta, by = "ID")
mean_expr_by_test$TEST = gsub("_CC", "", mean_expr_by_test$TEST)

png("Prental_Expression_Trajectories_CCa_CCo_min_100.png",
    res = 300,
    width = 15,
    height = 10,
    units = "in")

ggplot(mean_expr_by_test %>% na.omit() %>% filter(Period == "Prenatal"), 
       aes(x = Days, y = log(0.1 + Mean_RPKM), color = TEST, lty = TYPE, shape = TYPE)) + 
    geom_point() + 
    geom_smooth(method = "loess", formula = y ~ x, alpha = .1, span = 1) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom") + 
    geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(Region ~ .) +
    scale_color_manual(values = c("red",
                                  "blue",
                                  "green",
                                  "pink", 
                                  "yellow",
                                  "orange", 
                                  "black")) +
    scale_x_continuous(breaks = seq(0, 300, 15))

dev.off()

png("Postnatal_Expression_Trajectories_CCa_CCo_min_100.png",
    res = 300,
    width = 15,
    height = 10,
    units = "in")

ggplot(mean_expr_by_test %>% na.omit() %>% filter(Period == "Postnatal"), 
       aes(x = Days/365, y = log(0.1 + Mean_RPKM), color = TEST, lty = TYPE, shape = TYPE)) + 
    geom_point() + 
    geom_smooth(method = "loess", formula = y ~ x, alpha = .1, span = 1, se = FALSE) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom") + 
    geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(Region ~ .) +
    scale_color_manual(values = c("red",
                                  "blue",
                                  "green",
                                  "pink", 
                                  "yellow",
                                  "orange", 
                                  "black")) +
    scale_x_continuous(breaks = seq(0, 50, 5)) +
    xlab("Postnatal Years")

dev.off()

png("Pre_vsPostnatal_Expression_Trajectories_CCa_CCo_min_100.png",
    res = 300,
    width = 15,
    height = 10,
    units = "in")

ggplot(mean_expr_by_test %>% na.omit(), 
       aes(x = Days, y = log(0.1 + Mean_RPKM), color = TEST, lty = TYPE, shape = TYPE)) + 
    geom_point() + 
    geom_smooth(method = "loess", formula = y ~ x, alpha = .1, span = 1, se = FALSE) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom") + 
    geom_hline(yintercept = 0, lty = 2) +
    facet_grid(Region ~ TEST) +
    scale_color_manual(values = c("red",
                                  "blue",
                                  "green",
                                  "pink", 
                                  "yellow",
                                  "orange", 
                                  "black")) +
    scale_x_log10() +
    geom_vline(xintercept = 280, lty = 2)

dev.off()