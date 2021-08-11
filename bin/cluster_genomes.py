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
    default="average",
    help='Method used to perform linkage clustering (default is "average", which corresponds to UPGMA)'
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

# Function to read in a distmat matrix
def read_distmat(fp):
    logger.info("Reading .distmat file")

    # Make a list with the contents of each line
    output = []
    # Make a list with the index
    index = []

    # Iterate over the lines in the file
    for line_i, line in enumerate(open(fp, "r")):

        # If this is the first line
        if line_i == 0:

            # Skip it
            continue

        # For every other line
        else:

            # Remove the newline character at the end of the line
            line = line.rstrip("\n")

            # Remove every double space in the line
            while "  " in line:
                line = line.replace("  ", " ")

            # Split up the fields of the line
            fields = line.split(" ")
            
            # The index is the first position
            index.append(fields[0])

            # The rest of the fields are fields in the matrix
            output.append(fields[1:])

    # Format a DataFrame
    df = pd.DataFrame(
        output,
        index=index,
        columns=index
    ).applymap(
        float
    )

    # Return the DataFrame
    return df


# Read in the pairwise genome distances
logger.info(f"Reading from {args.dists}")
# If the distances were output by ClustalO
if args.dists.endswith(".distmat"):

    # Read the distances with a specific 'distmat' format
    dists = read_distmat(args.dists)

    # Get the name of the marker from the input file name
    msg = f"Input file does not conform to expected pattern: {args.dists}"
    assert args.dists.endswith(".markers.fasta.gz.distmat"), msg
    marker_name = args.dists.replace(".markers.fasta.gz.distmat", "")

    # Add the marker name to the output file name
    output_fp = f"{marker_name}.{args.ani_threshold}.hdf5"

# If the distances were output as a CSV
else:

    # Read in the distances in CSV format
    dists = pd.read_csv(
        args.dists,
        sep=",",
        index_col=0
    )

    # The output file will just contain the ANI in the filename
    output_fp = f"{args.ani_threshold}.hdf5"

logger.info(f"Read in {dists.shape[0]:,} rows and {dists.shape[1]:,} columns")
print(dists)

# Subset to genomes which have a valid distance measured
genome_list = [
    g for g in genome_list if g in dists.index.values
]
assert len(genome_list) > 0, "There are no genomes with alignments and the marker"

# Subset and order by the list of genomes with alignments
logger.info(f"Ordering distance matrix by the list of genomes")
dists = dists.reindex(
    index=genome_list,
    columns=genome_list,
).applymap(
    np.float16
)

# Make sure that the distances are symmetric
for genome_a in genome_list:
    dists.loc[genome_a, genome_a] = 0
    for genome_b in genome_list:
        if genome_a < genome_b:
            dists.loc[genome_a, genome_b] = dists.loc[genome_b, genome_a]




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
            n_genomes=group_df.genome.unique().shape[0],
            n_genes=group_df.sseqid.unique().shape[0]
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
with pd.HDFStore(output_fp, "w") as store:

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
