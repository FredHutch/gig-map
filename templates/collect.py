#!/usr/bin/env python3

import os
import pandas as pd


def read_tsv(fp):
    """
    The format of the mash dist outputs is typically:
    Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes

    In this workflow we have added the genome URI as the first column
    """

    return pd.read_csv(
        fp,
        sep="\\t",
        names = [
            "genome_uri",
            "ref",
            "query",
            "mash_distance",
            "p",
            "matching",
        ]
    )

# Read in all of the results
df = pd.concat([
    read_tsv(os.path.join('inputs', fp))
    for fp in os.listdir('inputs')
]).sort_values(
    by="mash_distance",
    ascending=True
)

# Write out the results
df.to_csv("${query_filename}.csv", sep=",", index=None)

# If the user wants to save any overlapping genomes
if "${params.save_overlapping}" != "false":

    # Get the number of minmers to use as a threshold
    try:
        min_val = int("${params.save_overlapping}")
    except:
        raise Exception("Could not parse as integer: ${params.save_overlapping}")

    # Filter the genomes by the minimum number of overlapping minmers
    overlapping = df.loc[
        df["matching"].apply(lambda s: s.split('/')[0]).apply(int) >= min_val
    ]

    # If there are any overlapping
    if overlapping.shape[0] > 0:

        # Write out just the genome URIs as a text file
        overlapping.reindex(columns=["genome_uri"]).to_csv("overlapping.txt", index=None, header=None)
