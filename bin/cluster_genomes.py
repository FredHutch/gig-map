#!/usr/bin/env python3
"""Cluster genomes by ANI."""

import argparse
import gzip
import logging
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

##################
# SET UP LOGGING #
##################

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [cluster_genomes.py] %(message)s'
)
logger = logging.getLogger('gig-map')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description="Cluster genomes by ANI"
)

# Add the arguments
parser.add_argument(
    '--alignments',
    type=str,
    required=True,
    help='Alignments of genes across genomes in CSV format'
)
parser.add_argument(
    '--dists',
    type=str,
    required=True,
    help='Pairwise ANI values for all genomes'
)
parser.add_argument(
    '--ani-threshold',
    type=int,
    required=True,
    help='Threshold used to cluster genomes (int: 0 - 100)'
)
parser.add_argument(
    '--method',
    type=str,
    default="ward",
    help='Method used to perform linkage clustering'
)

# Parse the arguments
args = parser.parse_args()

##############
# ALIGNMENTS #
##############

# Read in the alignments
logger.info(f"Reading from {args.alignments}")
alignments = pd.read_csv(args.alignments)

# Get a list of all genomes which have alignments
genome_list = list(alignments['genome'].unique())
logger.info(f"Read in a list of {len(genome_list):,} genomes which have alignments")


#############
# DISTANCES #
#############

# Read in the pairwise genome distances
logger.info(f"Reading from {args.dists}")
dists = pd.read_csv(
    args.dists,
    index_col=0
)
logger.info(f"Read in {dists.shape[0]:,} rows and {dists.shape[1]:,} columns")

# Subset and order by the list of genomes with alignments
logger.info(f"Ordering distance matrix by the list of genomes")
dists = dists.reindex(
    index=genome_list,
    columns=genome_list,
).applymap(
    np.float16
)



######################
# PERFORM CLUSTERING #
######################

# Start by performing linkage clustering
logger.info(f"Performing {args.method} linkage clustering on {dists.shape[0]} genomes")
L = linkage(squareform(dists.values), method=args.method)

# Get the flat clusters with this threshold
cluster_id_list = fcluster(
    L,
    1 - (args.ani_threshold / 100.),
    criterion="distance"
)

# Make a Series linking the genome to its group
genome_groups = pd.Series(
    cluster_id_list,
    index=dists.columns.values
)

# Reorder the indexes based on the size of each group
new_order = list(genome_groups.value_counts().index.values)
genome_groups = genome_groups.apply(
    new_order.index
)

logger.info("{} percent identity: Assigned {:,} genomes into {:,} groups".format(
    args.ani_threshold,
    genome_groups.shape[0],
    genome_groups.unique().shape[0]
))

# Make a summary of the genome groups
genome_group_df = pd.DataFrame(
    [
        dict(
            group=group_ix,
            n_genomes=group_df.genome_ix.unique().shape[0],
            n_genes=group_df.gene_ix.unique().shape[0]
        )
        for group_ix, group_df in alignments.groupby(
            alignments.genome.apply(
                genome_groups.get
            )
        )
    ]
# Sort the genome groups by the number of uniquely aligning genes
).sort_values(
    by="n_genes",
    ascending=False
).set_index(
    "group"
)

# Compute the average distances between genome groups
group_distance_df = dists.groupby(
    genome_groups
).mean(
).sort_index(
).T.groupby(
    genome_groups
).mean(
).sort_index(
)

################
# WRITE OUTPUT #
################

# Write out to an HDF5 file
with pd.HDFStore(f"{args.ani_threshold}.hdf5", "w") as store:

    # Write out the table of which genomes are assigned to which cluster
    logger.info("Writing out genome_groups")
    pd.DataFrame(
        [
            dict(
                genome=genome,
                group=group
            )
            for genome, group in genome_groups.items()
        ]
    ).to_hdf(store, "genome_groups")

    # Write out the table summarizing the contents of each cluster
    logger.info("Writing out genome_group_df")
    genome_group_df.to_hdf(store, "genome_group_df")

    # Write out the table of distances between each genome group
    logger.info("Writing out group_distance_df")
    group_distance_df.to_hdf(store, "group_distance_df")

logger.info("Done")
