#!/usr/bin/env python3
import pandas as pd
import gzip
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
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


def order_genes_large_scale(df):
    """Use PCA and t-SNE to set the order of a large number of genes."""

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
        n_components=1,
        n_jobs=n_cpus
    )

    # Group together genes which align to similar sets of genomes
    tsne_coords = tsne.fit_transform(
        pca_coords
    )

    # Get the gene order
    gene_order = pd.Series(
        tsne_coords[:,0],
        index=df.index.values
    ).sort_values(
    ).index.values

    return list(gene_order)


def order_genes_small_scale(df, method="ward", metric="euclidean"):
    """Use Ward linkage clustering to set the order of a small number of genes."""

    # Calculate the pairwise distances between genes
    dists = pdist(df.values, metric)

    # Perform linkage clustering
    Z = hierarchy.linkage(
        dists,
        method=method
    )

    # Order the tree so that adjacent leaves are more similar
    Z_ordered = hierarchy.optimal_leaf_ordering(
        Z,
        dists
    )

    # Get the ordered list of leaves
    leaves_list = hierarchy.leaves_list(Z_ordered)

    # Map that list of indexes back to the gene names
    gene_order = [
        df.index.values[i]
        for i in leaves_list
    ]

    return gene_order


# If there is a large number of genes
if df.shape[0] > max_n_genes_train_pca:

    # Use the PCA + t-SNE method used to perform ordination at large scale
    gene_order = order_genes_large_scale(df)

# If there is a smaller number of genes
else:

    # Use linkage clustering to perform ordination of genes
    gene_order = order_genes_small_scale(df)

# Write out to a file
with gzip.open(
    output_fp,
    "wt"
) as handle:
    handle.write("\n".join(gene_order))