#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(
    "unfiltered_gene_mapping.csv.gz",
    sep="\\t",
    header=None,
    names="${params.aln_fmt}".split(" ")
).rename(
    columns=dict(
        qseqid="query",
        sseqid="ref"
    )
).assign(
    coverage=lambda d: 100. * d['length'] / d['qlen'],
    score=lambda d: d['pident'] * d['coverage']
)

print(f"Read in {df.shape[0]:,} alignments")
print(f"Unique queries: {df['query'].unique().shape[0]:,}")
print(f"Unique references: {df['ref'].unique().shape[0]:,}")

print("Filtering by identity >= ${params.min_identity}")
df = df.query(
    "pident >= ${params.min_identity}"
)
print(f"Remaining alignments {df.shape[0]:,}")

print("Filtering by coverage >= ${params.min_coverage}")
df = df.query(
    "coverage >= ${params.min_coverage}"
)
print(f"Remaining alignments {df.shape[0]:,}")

print("Taking the top hit per query")
df = df.sort_values(
    by="score",
    ascending=False
).groupby(
    'query'
).head(
    1
).drop(
    columns=["score"]
)
print(f"Remaining alignments {df.shape[0]:,}")

df.to_csv("gene_mapping.csv.gz", sep=",", index=None)