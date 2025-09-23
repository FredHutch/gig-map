#!/usr/bin/env python3

import pandas as pd

# Read in the manifest
manifest = pd.read_csv("manifest.csv").set_index("specimen")
print(f"Read in {manifest.shape[0]:,} lines from the manifest")

# Read in the readcounts
read_alignments = pd.read_csv("read_alignments.csv.gz")
print(f"Read in {read_alignments.shape[0]:,} lines from the read_alignments")

# Filter the alignments to just those specimens in the manifest
read_alignments = read_alignments.loc[
    read_alignments["specimen"].isin(set(manifest.index.values))
]
print(f"After filtering to just those specimens in the manifest, {read_alignments.shape[0]:,} read alignments remain")
assert read_alignments.shape[0] > 0

# Make a list of all of the genes which have alignments
all_genes = read_alignments["id"].drop_duplicates().tolist()

# Break it up into shards of `shard_size` lengths
shard_size = int("${params.shard_size}")
gene_shards = [
    all_genes[i: min(i+shard_size, len(all_genes))]
    for i in range(0, len(all_genes), shard_size)
]

# If the final shard has only one member
if len(gene_shards[-1]) == 1 and len(gene_shards) > 1:
    gene_shards[-2].extend(gene_shards.pop())

# Make sure that all of the genes made it into a shard
assert len(all_genes) == sum(map(len, gene_shards)), (len(all_genes), sum(map(len, gene_shards)))

# Get the total number of reads per specimen
specimen_totals = read_alignments.groupby("specimen")["nreads"].sum()

# Write out each shard
for i, shard in enumerate(gene_shards):

    (
        read_alignments
        .loc[read_alignments["id"].isin(set(shard))]
        .assign(
            proportion=lambda d: d.apply(
                lambda r: r["nreads"] / specimen_totals[r['specimen']],
                axis=1
            )
        )
        .pivot_table(
            index="specimen",
            columns="id",
            values="proportion"
        )
        .reindex(index=manifest.index.values)
        .fillna(0)
        .to_csv(f"readcounts.{i}.csv.gz")
    )
