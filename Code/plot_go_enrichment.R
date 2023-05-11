#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(janitor)
})

list_of_files = list.files(path = getwd(), pattern = "*.gsa.out")
go_results = data.frame(GWAS = as.character(),
                        CATEGORY = as.character(), 
                        TERM = as.character(),
                        BETA = as.numeric(),
                        SE = as.numeric(),
                        P = as.numeric())

for (file in list_of_files) {
    go_out = fread(file, header = T)
    go_out = go_out %>% select(VARIABLE, FULL_NAME, BETA, SE, P) %>% 
        rename(CATEGORY = VARIABLE) %>%
        rename(TERM = FULL_NAME) %>%
        rename(P = P_ADJ) %>%
        mutate(GWAS = file)
    go_results = rbind(go_results, go_out)
}

# Housekeeping and cleaning

go_results$GWAS = gsub("^iPSYCH2015_EUR_", "", go_results$GWAS)
go_results$GWAS = gsub("_CC_C5_GO.gsa.out$", "", go_results$GWAS)
go_results$CATEGORY = gsub("_.*", "", go_results$CATEGORY)
go_results$FULL_NAME = make_clean_names(go_results$FULL_NAME)
go_results$P = go_results$P * length(list_of_files) # Bonferroni

# Choose unique enriched terms

enriched_go_terms = go_results %>% filter(P <= 0.05) %>% 
    select(TERM) %>% 
    unique()

go_enriched = inner_join(go_results, enriched_go_terms, by = c("TERM")) %>%
    mutate(Z = BETA/SE)

png("xDx2_GO_Enrichments.txt", 
    width = 12,
    height = 10,
    units = "in",
    res = 300)

ggplot(go_enriched, aes(x = TERM, y = GWAS, fill = Z)) +
    geom_tile() + 
    theme_classic() +
    geom_text(label = ifelse(P <= 0.05, "*", "")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
          axis.text.y = element_text(faace = "bold"),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    guides(fill = guide_colorbar(frame.colour = "black"))

dev.off()