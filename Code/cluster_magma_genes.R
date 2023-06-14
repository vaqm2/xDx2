#!/usr/bin/env Rscript

library(tidyverse)
library(factoextra)
library(gridExtra)
library(cluster)
library(dendextend)

set.seed(70783)
setwd("/Users/vapp0002/Documents/xDx_2/Magma/C_Magma_Genes/")

# Set number of clusters here

n_clusters = 2

# Initialize an empty data frame to collect magma associations

genes = data.frame(GENE = character(),
                   Z = numeric(),
                   P = numeric(),
                   TEST = character())
files = list.files(path = getwd(), pattern = "*.hgnc.out")

for (file in files) {
    print(paste("Processing", file, "....", sep = " "))
    association = read.table(file, header = T)
    test = gsub("^iPSYCH2015_EUR_", "", file)
    test = gsub("\\..*", "", test)
    association = association %>% 
        rename(Z = ZSTAT) %>% 
        mutate(TEST = test) %>% 
        select(GENE, Z, TEST, P)
    genes = rbind(genes, association)
}

genes_significant = genes %>% 
    group_by(TEST) %>%
    mutate(P_FDR = p.adjust(P, method = "fdr")) %>% 
    ungroup() %>%
    filter(P_FDR <= 0.05) %>%
    mutate(TYPE = str_count(TEST, "_"))

png("Genes_pass_FDR_by_test.png", 
    width = 8, 
    height = 8, 
    units = "in", 
    res = 300)

ggplot(genes_significant, aes(x = TEST)) + 
    geom_bar() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(. ~ TYPE, scales = "free", space = "free") +
    xlab("")

dev.off()

genes_significant = genes_significant %>%
    select(GENE) %>%
    unique()

genes_significant = inner_join(genes, genes_significant, by = "GENE") %>%
    mutate(TYPE = str_count(TEST, "_")) %>%
    select(-TYPE, -P) %>%
    reshape(direction = "wide", idvar = "GENE", timevar = "TEST")
colnames(genes_significant) = sub("Z.", "", colnames(genes_significant))
row.names(genes_significant) = genes_significant$GENE
genes_significant_scaled = genes_significant %>% select(-GENE) %>% scale()

# Computing distances between elements

gene_distance = get_dist(genes_significant_scaled, method = "pearson")

# Agglomerative hierarchical clustering

hc_agnes = agnes(gene_distance, method = "ward")

png("Agglomerative_hierarchical_clustering_dendrogram.png",
    width = 12,
    height = 12,
    res = 300,
    units = "in")

pltree(hc_agnes, cex = 0.6, hang = -1, main = "Agglomerative Clustering")

dev.off()

# Divisive hierarchical clustering

hc_diana = diana(gene_distance)

png("Divisive_hierarchical_clustering_dendrogram.png",
    width = 12,
    height = 12,
    res = 300,
    units = "in")

pltree(hc_diana, cex = 0.6, hang = -1, main = "Divisive Clustering")

dev.off()

# Clustering cut off evaluations

png("Elbow_method_evaluation.png",
    width = 8,
    height = 10,
    res = 300,
    units = "in")

fviz_nbclust(genes_significant_scaled, 
             FUN = hcut, 
             method = "wss", 
             k.max = 100)

dev.off()

png("Silhouette_method_evaluation.png",
    width = 10,
    height = 10,
    res = 300,
    units = "in")

fviz_nbclust(genes_significant_scaled, 
             FUN = hcut, 
             method = "silhouette", 
             k.max = 100)

dev.off()

gap_stat = clusGap(genes_significant_scaled, 
                   FUN = hcut, 
                   nstart = 25, 
                   K.max = 100, 
                   B = 50)

png("Gap_stat_evluation.png",
    width = 10,
    height = 10,
    res = 300,
    units = "in")

fviz_gap_stat(gap_stat)

dev.off()

groups_agnes = cutree(hc_agnes, k = n_clusters) %>% as.data.frame()
colnames(groups_agnes) = c("GROUP")
groups_agnes$GENE = row.names(groups_agnes)
genes = inner_join(genes, groups_agnes, by = c("GENE"))

mean_z = genes %>% 
    select(-GENE) %>% 
    group_by(GROUP, TEST) %>% 
    mutate(MEAN_Z = mean(Z)) %>% 
    select(-Z) %>% 
    unique() %>% 
    ungroup() %>% 
    mutate(TYPE = str_count(TEST, "_"))

mean_z$TYPE = gsub("0", "[a] Case vs. Cohort", mean_z$TYPE)
mean_z$TYPE = gsub("1", "[b] Case vs. Other Cases", mean_z$TYPE)
mean_z$TYPE = gsub("2", "[c] Case vs. Case Pairwise", mean_z$TYPE)
mean_z$TEST = gsub("_CC", "", mean_z$TEST)
mean_z$TEST = gsub("_", " vs. ", mean_z$TEST)

png(paste0("Mean_Z_by_cluster_group", n_clusters, ".png"),
    width = 12,
    height = 8,
    res = 300,
    units = "in")

ggplot(mean_z, aes(y = factor(GROUP), x = TEST, fill = MEAN_Z)) + 
    geom_tile(color = "black") + 
    geom_text(aes(label = round(MEAN_Z, 2)), size = 3) + 
    theme_classic() + 
    facet_grid(. ~ TYPE, scales = "free", space = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    ylab("CLUSTER") + 
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"), 
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom")

dev.off()

png(paste0("Num_genes_by_cluster_group", "_", num_clusters, ".png"),
    res = 300,
    units = "in",
    width = 10,
    height = 10)

ggplot(mean_z, aes(x = factor(GROUP), fill = factor(GROUP))) + 
    geom_bar() + 
    theme_classic() +
    guides(fill = "none") +
    xlab("CLUSTER")

dev.off()