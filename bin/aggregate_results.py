#!/usr/bin/env python3
"""Script to aggregate all results of the gig-map processing for rapid visualization."""

import argparse
from direct_redis import DirectRedis
import gzip
import logging
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


#############
# DISTANCES #
#############

# Read in the pairwise genome distances
logger.info(f"Reading from {args.dists}")
dists = pd.read_csv(
    args.dists,
    index_col=0
)


#####################
# t-SNE COORDINATES #
#####################

# Read in the t-SNE coordinates per-gene
logger.info(f"Reading from {args.tnse_coords}")
tsne_coords = pd.read_csv(
    args.tnse_coords,
    index_col=0
)


###############
# GENOME LIST #
###############

# Get a list of all genomes which have alignments
genome_list = list(alignments['genome'].unique())


#############
# GENE LIST #
#############

# Read the list of all genes, ordered by similarity of alignment
gene_list = [
    line.decode().rstrip("\n")
    for line in gzip.open(args.gene_order, 'r')
]

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

    # Save the alignment information in a single table
    logger.info("Saving alignments to redis")
    r.set(
        # The key at which the values may be accessed
        "alignments",
        # The values which will be accessed at the key
        alignments
    )

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

    # Save the table of distances
    logger.info("Saving distances to redis")
    r.set(
        # The key at which the values may be accessed
        "distances",
        # The values which will be accessed at the key
        dists
    )
    
    # Save the table of t-SNE coordinates
    logger.info("Saving tsne to redis")
    r.set(
        # The key at which the values may be accessed
        "tsne",
        # The values which will be accessed at the key
        tsne_coords
    )
