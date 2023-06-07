#!/usr/bin/env python

import os
import re
import sys
import pandas as pd
import glob
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from matplotlib import pyplot as plt
import seaborn as sn

data_dir = sys.argv[1]
magma_gene_files = glob.glob(data_dir + '/*.hgnc.out')
magma_gene_dfs = []

for file in magma_gene_files:
    print(f'Reading file: {file}')
    trait = os.path.basename(file)
    trait = re.sub("^iPSYCH2015_EUR_", "", trait)
    trait = re.sub("\..*$", "", trait)
    file_df = pd.read_csv(file,
                          sep = '\s+',
                          usecols = ['GENE', 'ZSTAT'], 
                          index_col = ['GENE']).rename(columns = {'ZSTAT' : trait}, 
                                                       inplace = False)
    magma_gene_dfs.append(file_df)

assoc_df = pd.concat(magma_gene_dfs, axis = 1)

cluster_complete = linkage(assoc_df, 
                           method = "complete", 
                           metric = "euclidean")

# dendrogram(cluster_complete)
# plt.savefig('cluster_complete_genes.png')
clusters = cut_tree(cluster_complete, 7)
assoc_df = assoc_df.assign(CLUSTER = clusters)
cluster_means = assoc_df.groupby(['CLUSTER']).mean()
cluster_means_df = pd.DataFrame(cluster_means)

fig, ax = plt.subplots(figsize=(14, 6))

sn.heatmap(data = cluster_means_df,
           cmap = 'coolwarm',
           annot = True,
           fmt = '.2g',
           cbar = True,
           square = True,
           linewidths = 0.5,
           linecolor = 'white')

plt.show()