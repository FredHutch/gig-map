#!/usr/bin/env python3

import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import logging
import sys
logger = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Get the threshold used for filtering
logger.info("Filtering ANI values below (filter_ani_threshold): ${params.filter_ani_threshold}")
filter_ani_threshold = float("${params.filter_ani_threshold}")

# Read in the table of distances between genomes
dists = pd.read_csv("distances.csv.gz", index_col=0)
logger.info("Read in the distances table: %s", dists.shape)

# Run average linkage clustering
linkage_matrix = linkage(squareform(dists), method='average')

# Get the flat clusters at that threshold
clusters = fcluster(linkage_matrix, filter_ani_threshold, criterion='distance')
logger.info("Got %d clusters", len(set(clusters)))

# Find the largest cluster
clusters_vc = pd.Series(clusters).value_counts()
largest_cluster = clusters_vc.idxmax()
logger.info("Largest cluster is %d (%d members)", largest_cluster, clusters_vc.max())

# Write out a list of the IDs of the  genomes in the largest cluster
with open("filtered_genomes_ids.txt", "w") as ofh:
    ofh.write("\\n".join([
        genome_id
        for genome_id, cluster in zip(dists.index, clusters)
        if cluster == largest_cluster
    ]))
logger.info("Wrote the list of genomes in the largest cluster to filtered_genomes_ids.txt")
