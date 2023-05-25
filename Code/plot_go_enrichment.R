#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(janitor)
    library(scales)
    library(stringr)
})

list_of_files = list.files(path = getwd(), pattern = "*_sets.gsa.out")
go_results = data.frame(GWAS = as.character(),
                        CATEGORY = as.character(), 
                        TERM = as.character(),
                        BETA = as.numeric(),
                        SE = as.numeric(),
                        P = as.numeric())

for (file in list_of_files) {
    go_out = fread(file, skip = 2, header = T)
    go_out = go_out %>% 
        filter(NGENES >= 10 & NGENES <= 500) %>%
        select(VARIABLE, FULL_NAME, BETA, SE, P) %>% 
        rename(CATEGORY = VARIABLE) %>%
        rename(TERM = FULL_NAME) %>%
        mutate(GWAS = file)
    go_results = rbind(go_results, go_out)
}

# Housekeeping and cleaning

go_results$GWAS = gsub("^iPSYCH2015_EUR_", "", go_results$GWAS)
go_results$GWAS = gsub("_C5_GO.gsa.out$", "", go_results$GWAS)
go_results$CATEGORY = gsub("_.*", "", go_results$CATEGORY)
go_results$CATEGORY = gsub("^GO", "", go_results$CATEGORY)
go_results$P = go_results$P * length(list_of_files) # Bonferroni
go_results = go_results %>% mutate(TYPE = str_count(GWAS, "_"))

# Choose unique enriched terms

enriched_go_terms = go_results %>% 
    filter(P <= 0.05 & TYPE < 2) %>% 
    select(TERM) %>% 
    unique()
enriched_go_terms_pairwise = go_results %>% 
    filter(P <= 0.05 & TYPE == 2) %>% 
    select(TERM) %>% 
    unique()
go_enriched = inner_join(go_results, enriched_go_terms, by = c("TERM")) %>%
    filter(TYPE < 2) %>%
    mutate(Z = BETA/SE)
go_enriched_pairwise = inner_join(go_results, enriched_go_terms_pairwise, 
                                  by = c("TERM")) %>%
    filter(TYPE == 2) %>%
    mutate(Z = BETA/SE)
go_enriched = rbind(go_enriched, go_enriched_pairwise)
    

# More housekeeping and cleaning

go_enriched$TERM = gsub("^.*?_", "", go_enriched$TERM)
go_enriched$TERM = make_clean_names(go_enriched$TERM, allow_dupes = TRUE)
go_enriched$TERM = gsub("_", " ", go_enriched$TERM)
go_enriched$TYPE = gsub("0", "Case vs Cohort", go_enriched$TYPE)
go_enriched$TYPE = gsub("1", "Case vs Other Cases", go_enriched$TYPE)
go_enriched$TYPE = gsub("2", "Case vs Case Pairwise", go_enriched$TYPE)
go_enriched$GWAS = gsub("_CC", "", go_enriched$GWAS)
go_enriched$GWAS = gsub("_", " vs ", go_enriched$GWAS)

png("xDx2_Konrad_Enrichments.png", 
    width = 12,
    height = 15,
    units = "in",
    res = 300)

ggplot(go_enriched %>% filter(TYPE != "Case vs Case Pairwise"), 
       aes(y = TERM, x = GWAS, fill = Z)) +
    geom_tile(color = "black") + 
    theme_classic() +
    geom_text(aes(label = ifelse(P <= 0.05, "*", ""))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    guides(fill = guide_colorbar(frame.colour = "black")) +
    scale_y_discrete(labels = label_wrap(60), guide = guide_axis(n.dodge = 2)) +
    facet_grid(. ~ TYPE, scales = "free", space = "free") +
    ylab("") + xlab("")

dev.off()

png("xDx2_Konrad_Enrichments_pairwise.png", 
    width = 12,
    height = 15,
    units = "in",
    res = 300)

ggplot(go_enriched %>% filter(TYPE == "Case vs Case Pairwise"), 
       aes(y = TERM, x = GWAS, fill = Z)) +
    geom_tile(color = "black") + 
    theme_classic() +
    geom_text(aes(label = ifelse(P <= 0.05, "*", ""))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    guides(fill = guide_colorbar(frame.colour = "black")) +
    scale_y_discrete(labels = label_wrap(60), guide = guide_axis(n.dodge = 2)) +
    facet_grid(. ~ TYPE, scales = "free", space = "free") +
    ylab("") + xlab("")

dev.off()