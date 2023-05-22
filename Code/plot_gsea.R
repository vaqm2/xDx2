#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(janitor)
    library(stringr)
    library(reshape2)
})

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

for (file in gsea_files) {
    gwas = gsub("_GSEA_Results.txt$", "", file)
    gwas = gsub("_kegg|_reactome|_wikipathways|_bp|_mf|_cc", "", gwas)
    gsea = fread(file, header = T)
    gsea = gsea %>% 
        select(pathway, pval, padj, ES, NES) %>% 
        mutate(study = gwas)
    gsea_results = rbind(gsea_results, gsea)
}

gsea_results = gsea_results %>% mutate(Type = str_count(study, "_"))
gsea_results$Type = gsub("0", "Case vs Cohort", gsea_results$Type)
gsea_results$Type = gsub("1", "Case vs Other Cases", gsea_results$Type)
gsea_results$Type = gsub("2", "Case vs Case Pairwise", gsea_results$Type)
gsea_results$study = gsub("_CC", "", gsea_results$study)
gsea_results$study = gsub("_", " vs ", gsea_results$study)
gsea_results_xdx = gsea_results %>% filter(study == "xDx")
gsea_results_xdx$Type = "xDx"
gsea_results = gsea_results %>% filter(study != "xDx")
gsea_results = rbind(gsea_results, gsea_results_xdx)

for (source in c("WP")) {
    data = gsea_results %>% filter(str_starts(pathway, source))
    data = inner_join(data, enriched_pathways, by = c("pathway"))
    data$pathway = make_clean_names(data$pathway, allow_dupes = TRUE)
    data$pathway = gsub("^gobp_|^gomf_|^gocc_|^kegg_|^wp_|^reactome_", "", 
                        data$pathway)
    data$pathway = gsub("_", " ", data$pathway)
    print(
        ggplot(data, aes(y = pathway, x = study, fill = NES)) + 
            geom_tile(color = "white") +
            theme_classic() + 
            geom_text(aes(label = ifelse(padj <= 0.05, "*", ""))) +
            scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                  axis.text.y = element_text(face = "bold")) + 
            xlab("") + 
            ylab("") + 
            facet_grid(. ~ Type, scales = "free", space = "free") +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 40)))
}