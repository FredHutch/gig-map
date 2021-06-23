#!/usr/bin/env python3
"""Script to aggregate all results of the gig-map processing for rapid visualization."""

import gzip
import pandas as pd
import argparse
import logging
import pickle
pickle.HIGHEST_PROTOCOL = 4

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
    '--output',
    type=str,
    required=True,
    help='Location to write HDF5 store with aggregated data'
)

# Parse the arguments
args = parser.parse_args()

# Make sure that the output path has the expected extension
assert args.output.endswith(".hdf5"), "--output must end with .hdf5"

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

def remove_genome_file_ext(fp):
    for ext in ['.gz', '.fna', '.fasta', '.fa']:
        if fp.endswith(ext):
            fp = fp[:-len(ext)]
    return fp

# Remove the file endings from the genome names
alignments = alignments.apply(
    lambda c: c.apply(remove_genome_file_ext) if c.name == "genome" else c
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
    sep="\t",
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
genome_list = alignments['genome'].unique()


#############
# GENE LIST #
#############

# Read the list of all genes, ordered by similarity of alignment
gene_list = [
    line.decode().rstrip("\n")
    for line in gzip.open(args.gene_order, 'r')
]


################
# WRITE OUTPUT #
################

# Open a connection to the output HDF store
with pd.HDFStore(args.output, 'w') as store:
    # Save the alignment information in three tables
    # All three tables will have the same index and columns
    # The index will be the list of genomes
    # The columns will be the genes, ordered by `--gene-order`
    # /alignments/pident
    #       Every value is the maximum percent identity of alignment
    # /alignments/coverage
    #       Every value is the maximum coverage of alignment
    # /alignments/description
    #       Every value is a text string summarizing the location of all alignments

    # Iterate over each column in the deduplicated table
    for colname in alignments.columns.values:

        alignments.pivot(
            index="genome",
            columns="sseqid",
            values=colname
        ).reindex(
            index=genome_list,
            columns=gene_list
        ).to_hdf(
            store,
            f"/alignments/{colname}",
            complevel=9,
            format="fixed",
        )

    # Save the table of distances
    dists.to_hdf(
        store,
        "/distances",
        complevel=9,
        format="fixed"
    )
    
    # Save the table of t-SNE coordinates
    tsne_coords.to_hdf(
        store,
        "/tsne",
        complevel=9,
        format="fixed"
    )