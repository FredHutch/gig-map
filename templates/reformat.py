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
            "ref",
            "query",
            "mash_distance",
            "p",
            "matching",
        ]
    )

# Read in all of the results
for fp in os.listdir("."):
    if fp.endswith('.tsv'):
        df = read_tsv(fp)

# Sort by distance
df.sort_values(
    by="mash_distance",
    ascending=True,
    inplace=True
)

# Remove the .orfs.fasta from the genomes (if any)
df = df.assign(
    ref = df["ref"].apply(
        lambda s: s[:-len(".orfs.fasta") if s.endswith(".orfs.fasta") else s]
    )
)

# Save a table for each query
for query_name, query_df in df.groupby('query'):
    query_df.to_csv(f"{query_name}.dists.csv", index=None)
