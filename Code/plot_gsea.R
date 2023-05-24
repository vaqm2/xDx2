#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(janitor)
    library(stringr)
    library(ggrepel)
})

# Working directory and glob fetch files by pattern 

setwd("/Users/vapp0002/Documents/xDx_2/Magma/GSEA/")
gsea_files = list.files(pattern = "*_GSEA_Results.txt")
gsea_files_concise = list.files(pattern = "*_GSEA_Concise_Results.txt")

enriched_pathways = data.frame(pathway = character())
gsea_results = data.frame(study = character(), 
                          pathway = character(),
                          pval = integer(),
                          padj = integer(),
                          ES = integer(),
                          NES = integer(),
                          stringsAsFactors = FALSE)

# Find all unique concise pathway names

for (file in gsea_files_concise) {
    gsea = fread(file, header = T)
    gsea_pathways = gsea %>% 
        filter(padj <= 0.05) %>% 
        select(pathway)
    enriched_pathways = rbind(enriched_pathways, gsea_pathways)
}

enriched_pathways = enriched_pathways %>% 
    select(pathway) %>% 
    unique()

# Read in all the pathway test results

for (file in gsea_files) {
    gwas = gsub("_GSEA_Results.txt$", "", file)
    gwas = gsub("_kegg|_reactome|_wikipathways|_bp|_mf|_cc", "", gwas)
    gsea = fread(file, header = T)
    gsea = gsea %>% 
        select(pathway, pval, padj, ES, NES) %>% 
        mutate(study = gwas)
    gsea_results = rbind(gsea_results, gsea)
}

# Subset to concise pathways, label and clean studies

gsea_results = inner_join(gsea_results, enriched_pathways, by = c("pathway")) %>%
    mutate(Resource = str_replace(pathway, "_.*", ""))
gsea_results = gsea_results %>% mutate(Type = str_count(study, "_"))
gsea_results$study = gsub("_CC", "", gsea_results$study)
gsea_results$pathway = make_clean_names(gsea_results$pathway, allow_dupes = TRUE)
gsea_results$pathway = gsub("^gobp_|^gomf_|^gocc_|^kegg_|^wp_|^reactome_", "", 
                        gsea_results$pathway)

# Separate cross disorder/pleiotropy results into a data frame

gsea_results_xdx = gsea_results %>% 
    filter(study == "xDx") %>% 
    select(pathway, padj, NES, Resource) %>%
    rename(padj_xdx = padj) %>% 
    rename(NES_xdx = NES)

# Separate case-control pathway enrichments into a data frame

gsea_results_cco = gsea_results %>% 
    filter(Type == 0 & study != "xDx") %>% 
    select(pathway, padj, NES, study, Resource) %>%
    rename(padj_cco = padj) %>% 
    rename(NES_cco = NES)

# Separate case-case pathway enrichments into a data frame

gsea_results_cca = gsea_results %>% 
    filter(Type == 1) %>%
    select(pathway, padj, NES, study, Resource) %>%
    rename(padj_cca = padj) %>% 
    rename(NES_cca = NES)

# and finally the pairwise results

gsea_results_pairwise = gsea_results %>% 
    filter(Type != 1 & study != "xDx") %>% 
    select(pathway, padj, NES, study, Resource)

# INTERSECT results for plotting

# Case-Case vs Cross Disorder

gsea_results_cca_xdx = inner_join(gsea_results_cca, gsea_results_xdx, 
                                  by = c("pathway", "Resource")) %>%
    mutate(Enriched = case_when(padj_cca <= 0.05 & padj_xdx <= 0.05 ~ "Both",
                                padj_cca <= 0.05 & padj_xdx > 0.05 ~ "Specific",
                                padj_cca > 0.05 & padj_xdx <= 0.05 ~ "Pleiotropic",
                                .default = "Neither"))

ggplot(gsea_results_cca_xdx, aes(x = NES_xdx, y = NES_cca, color = Enriched)) + 
    geom_point(pch = 21) + 
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope = 1, lty = 2, color = "black") + 
    facet_grid(Resource ~ study, scales = "free", space = "free") + 
    scale_color_manual(values = c("green", "gray", "blue", "red")) +     
    xlab("Normalized Enrichment Score Cross Disorder") +
    ylab("Normalized Enrichment Score Case vs Other Cases") +     
    geom_label_repel(aes(label = ifelse(Enriched == "Specific", pathway, "")),
                    size = 3,
                    max.overlaps = Inf)

# Case-Case vs Case-Control

gsea_results_cca_cco = inner_join(gsea_results_cca, gsea_results_cco, 
                                  by = c("pathway", "study", "Resource")) %>%
    mutate(Enriched = case_when(padj_cca <= 0.05 & padj_cco <= 0.05 ~ "Both",
                                padj_cca <= 0.05 & padj_cco > 0.05 ~ "vs Cohort",
                                padj_cca > 0.05 & padj_cco <= 0.05 ~ "vs Other Cases",
                                .default = "Neither"))

ggplot(gsea_results_cca_cco, aes(x = NES_cco, y = NES_cca, color = Enriched)) + 
    geom_point(pch = 21) + 
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope = 1, lty = 2, color = "black") + 
    facet_wrap(study ~ .) + 
    scale_color_manual(values = c("green", "gray", "blue", "red")) +
    xlab("Normalized Enrichment Score Case vs Cohort") +
    ylab("Normalized Enrichment Score Case vs Other Cases") +
    geom_text_repel(aes(label = ifelse(Enriched == "vs Other Cases", pathway, "")), 
                    size = 3,
                    max.overlaps = Inf)