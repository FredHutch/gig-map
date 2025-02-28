#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram

# Read in all of the distances
df = pd.concat(
    [
        pd.read_csv(
            os.path.join('inputs', fp),
            sep="\\t",
            header=None,
            names=['query', 'ref', 'dist', 'n', 'ratio']
        ).reindex(
            columns=['query', 'ref', 'dist']
        )
        for fp in os.listdir('inputs')
    ]
).pivot_table(
    index='query',
    columns='ref',
    values='dist'
)

# Save the distances to CSV
df.to_csv("distances.csv.gz")

# Check to see if the filter_ani_threshold parameter is set
try:
    filter_ani_threshold = float("${params.filter_ani_threshold}")
except: # noqa
    filter_ani_threshold = None

# Build a linkage matrix using this distance matrix
linkage_matrix = linkage(squareform(df), method='average')

# Color the dendrogram based on the filter_ani_threshold
dendrogram(
    linkage_matrix,
    color_threshold=filter_ani_threshold,
    labels=df.index,
    orientation='left'
)
# Adjust the height of the plot to fit the labels
plt.gcf().set_size_inches(8, 0.15 * len(df.index))
plt.xlabel("ANI Distance")
try:
    plt.tight_layout()
except:
    pass
try:
    plt.savefig('ani_dendrogram.png')
    plt.savefig('ani_dendrogram.pdf')
except:
    pass
