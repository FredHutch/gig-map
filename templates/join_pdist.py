#!/usr/bin/env python3

from collections import defaultdict
import gzip
import json
import os
import pandas as pd

# Format the output as 
# {
#     min_name: {
#         max_name: pident
#     }
# }


def read_alignment(fp):
    df = pd.read_csv(
        fp,
        sep="\\t",
        header=None,
        names="${params.aln_fmt}".split(" ")
    ).rename(
        columns=dict(
            qseqid="query",
            sseqid="ref"
        )
    ).assign(
        coverage=lambda d: 100. * d['length'] / d['qlen']
    )

    print(f"Read in {df.shape[0]:,} alignments from {fp}")
    
    print("Removing self alignments")
    df = df.query("query != ref")
    print(f"Remaining alignments {df.shape[0]:,}")

    print("Filtering by coverage >= ${params.min_coverage}")
    df = df.query(
        "coverage >= ${params.min_coverage}"
    )
    print(f"Remaining alignments {df.shape[0]:,}")

    for _, r in df.iterrows():
        min_name = min(r['query'], r['ref'])
        max_name = max(r['query'], r['ref'])
        yield min_name, max_name, r['pident']

output = defaultdict(dict)

# Iterate over each of the input files
for input_fp in os.listdir("."):
    if not input_fp.endswith(".gene_mapping.csv.gz"):
        continue
    for min_name, max_name, pident in read_alignment(input_fp):
        output[min_name][max_name] = pident

# Write out the data
print("Writing out data")
with gzip.open("gene_pdist.json.gz", "wt") as handle:
    json.dump(output, handle)
print("Done")
