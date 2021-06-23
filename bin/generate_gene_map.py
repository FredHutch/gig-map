#!/usr/bin/env python3
import pandas as pd
import gzip
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import sys

# Get the user-supplied arguments
input_fp = sys.argv[1]
max_n_genes_train_pca = int(sys.argv[2])
max_pcs_tsne = int(sys.argv[3])
output_fp = sys.argv[4]
n_cpus = int(sys.argv[5])

# Read in all of the alignment information
df = pd.read_csv(
    input_fp
# Pivot to wide format
).pivot_table(
    columns="genome",
    index="sseqid",
    values="pident"
).fillna(
    0
)

# Train a PCA projection using a subset of the data

# Initialize the PCA object
pca = PCA(
    n_components=min(
        df.shape[1],
        max_pcs_tsne
    )
)

# Fit with a subset of the genes
pca.fit(
    df.sample(
        min(
            df.shape[0],
            max_n_genes_train_pca
        )
    )
)

# Get the PCA coordinates for the full set of genes
pca_coords = pca.transform(
    df.values
)

# Initialize the TSNE object
tsne = TSNE(
    n_components=2,
    n_jobs=n_cpus
)

# Group together genes which align to similar sets of genomes
tsne_coords = tsne.fit_transform(
    pca_coords
)

# Format a DataFrame with the results
gene_map_df = pd.DataFrame(
    tsne_coords,
    index=df.index.values,
    columns=[
        "t-SNE 1",
        "t-SNE 2",
    ]
)

# Write out to a file
gene_map_df.to_csv(
    output_fp,
    compression="gzip"
)
