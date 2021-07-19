#!/usr/bin/env python3
"""Script to aggregate all results of the gig-map processing for rapid visualization."""

import argparse
from direct_redis import DirectRedis
import gzip
import logging
import os
import numpy as np
import pandas as pd

##################
# SET UP LOGGING #
##################

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [aggregate_results.py] %(message)s'
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
    description="Aggregate all results of the gig-map processing for rapid visualization"
)

# Add the arguments
parser.add_argument(
    '--alignments',
    type=str,
    required=True,
    help='Alignments of genes across genomes in CSV format'
)
parser.add_argument(
    '--gene-order',
    type=str,
    required=True,
    help='Ordering of genes by presence across genomes'
)
parser.add_argument(
    '--dists',
    type=str,
    required=True,
    help='Pairwise ANI values for all genomes'
)
parser.add_argument(
    '--tnse-coords',
    type=str,
    required=True,
    help='t-SNE coordinates for all genes in two dimensions'
)
parser.add_argument(
    '--host',
    type=str,
    default="localhost",
    help='Redis host used for writing output'
)
parser.add_argument(
    '--port',
    type=int,
    default=6379,
    help='Redis port used for writing output'
)
parser.add_argument(
    '--dists-n-rows',
    type=int,
    default=1000,
    help='Number of rows to use for each chunk of distances'
)
parser.add_argument(
    '--alignments-n-rows',
    type=int,
    default=1000000,
    help='Number of rows to use for each chunk of alignments'
)

# Parse the arguments
args = parser.parse_args()

##############
# ALIGNMENTS #
##############

# Read in the alignments
logger.info(f"Reading from {args.alignments}")
alignments = pd.read_csv(args.alignments)

# Calculate the alignment coverage of each gene
alignments = alignments.assign(
    coverage = alignments.apply(
        lambda r: 100 * (r['send'] - r['sstart'] + 1) / r['slen'],
        axis=1
    )
)


#######################
# REFORMAT ALIGNMENTS #
#######################

# Set up a function to format a text string describing >=1 alignments
def format_description(d):
    
    return "\n".join([
        f"{r['qseqid']}: {r['qstart']:,} - {r['qend']:,}; {r['pident']}% identity / {r['coverage']}% coverage"
        for _, r in d.iterrows()
    ])

# Condense the alignment table to just have a single row per gene/genome
alignments = alignments.groupby(
    ["sseqid", "genome"]
).apply(
    lambda d: pd.Series(
        dict(
            pident=d['pident'].max(),
            coverage=d['coverage'].max(),
            description=format_description(d)
        )
    )
).reset_index()


#####################
# t-SNE COORDINATES #
#####################

# Read in the t-SNE coordinates per-gene
logger.info(f"Reading from {args.tnse_coords}")
tsne_coords = pd.read_csv(
    args.tnse_coords,
    index_col=0
)
logger.info(f"Read in {tsne_coords.shape[0]:,} rows and {tsne_coords.shape[1]:,} columns")


###############
# GENOME LIST #
###############

# Get a list of all genomes which have alignments
genome_list = list(alignments['genome'].unique())
logger.info(f"Read in a list of {len(genome_list):,} genomes")


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


###################
# GENOME CLUSTERS #
###################

# Read in each of the groupings which result from ANI-based genome clustering
genome_clustering = dict()

# Each set of results is available as genome_clusters/{ani_threshold}.hdf5
for fp in os.listdir("genome_clusters"):

    # Only read in HDF5 files
    if not fp.endswith(".hdf5"):
        continue

    # Parse the threshold from the file name
    ani_threshold = fp[:-5]
    logger.info(f"Parsing results of clustering at {ani_threshold} percent identity")
    
    # The threshold must be an integer (indicating a percentage)
    ani_threshold = int(ani_threshold)

    # Open the HDF5 file for reading
    fp = f"genome_clusters/{fp}"
    logger.info(f"Opening {fp}")
    with pd.HDFStore(fp, "r") as store:

        # Iterate over a set of keys
        for hdf_key in ["genome_groups", "genome_group_df", "group_distance_df"]:

            logger.info(f"Reading key {hdf_key}")

            # Save the DataFrame to `genome_clustering`
            genome_clustering[
                f"{hdf_key} {ani_threshold}"
            ] = pd.read_hdf(
                store,
                hdf_key
            )

logger.info("Done reading genome clusters")


#############
# GENE LIST #
#############

# Read the list of all genes, ordered by similarity of alignment
gene_list = [
    line.decode().rstrip("\n")
    for line in gzip.open(args.gene_order, 'r')
]
logger.info(f"Read in a list of {len(gene_list):,} genes")

# Add the index position of each gene and genome to the table
alignments = alignments.assign(
    gene_ix = alignments.sseqid.apply(
        gene_list.index
    ),
    genome_ix = alignments.genome.apply(
        genome_list.index
    )
).drop(
    columns=["sseqid", "genome"]
)


################
# WRITE OUTPUT #
################

# Connect to redis
logger.info(f"Connecting to redis at {args.host}:{args.port}")
with DirectRedis(host=args.host, port=args.port) as r:

    # Save the alignment information in a set of tables
    logger.info(f"Saving {alignments.shape[0]:,} alignments to redis")
    alignments_keys = []

    # Iterate over chunks of the DataFrame
    for chunk_ix, chunk_start in enumerate(np.arange(0, alignments.shape[0], args.alignments_n_rows)):

        # Set the ending row index of the chunk
        chunk_end = min(chunk_start + args.alignments_n_rows, alignments.shape[0])

        # Set the key used for this chunk in redis
        chunk_key = f"alignments_{chunk_ix}"

        # Write the chunk of the DataFrame
        r.set(
            # The key at which the values may be accessed
            chunk_key,
            # The values which will be accessed at the key
            alignments.loc[
                chunk_start: chunk_end
            ]
        )

        # Add the key to the list
        alignments_keys.append(chunk_key)

    # Save the list of chunk keys to redis
    r.set("alignments_keys", alignments_keys)

    # Save the mapping of gene_ix to a name
    logger.info("Saving gene_ix to redis")
    r.set(
        # The key at which the values may be accessed
        "gene_ix",
        # The values which will be accessed at the key
        gene_list
    )

    # Save the mapping of genome_ix to a name
    logger.info("Saving genome_ix to redis")
    r.set(
        # The key at which the values may be accessed
        "genome_ix",
        # The values which will be accessed at the key
        genome_list
    )

    # Save the table of distances in chunks of `dists_n_rows` each
    logger.info("Saving distances to redis")

    # Keep track of the keys used to store chunks in the database
    dists_keys = []

    # Iterate until all of the distances have been written
    while dists.shape[0] > 0:

        # Write a chunk of distances
        r.set(
            # Key the chunk by the index
            f"distances_{len(dists_keys)}",
            # Write the first `dists_n_rows` rows to redis
            dists.iloc[:min(dists.shape[0], args.dists_n_rows)]
        )

        # Keep track of the key that was used
        dists_keys.append(f"distances_{len(dists_keys)}")

        logger.info(f"Wrote {len(dists_keys):,} chunks of distances")

        # If the complete set of distances has been written
        if dists.shape[0] <= args.dists_n_rows:

            # Stop iterating
            break

        # Otherwise
        else:

            # Remove the written chunk from the overall dists table
            dists = dists.iloc[args.dists_n_rows:]

    logger.info(f"Done writing distances")

    # Store the list of keys which were used for the chunks
    r.set(
        "distances_keys",
        dists_keys
    )

    # Write out all of the genome clustering information
    for k, v in genome_clustering.items():
        logger.info(f"Writing {k} to redis")
        r.set(k, v)
    
    # Save the table of t-SNE coordinates
    logger.info("Saving tsne to redis")
    r.set(
        # The key at which the values may be accessed
        "tsne",
        # The values which will be accessed at the key
        tsne_coords
    )

    logger.info("Done writing to redis")

logger.info("Closed connection to redis")
