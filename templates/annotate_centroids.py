#!/usr/bin/env python3

import gzip
import pandas as pd
import os

# Read all of the annotations from the input files
gene_annot = dict()

for fp in os.listdir("input_annotations"):
    print(f"Reading in annotations from {fp}")

    for _, r in pd.read_csv(
        os.path.join("input_annotations", fp)
    ).iterrows():
        if r.get("gene_id") is not None and r.get("combined_name") is not None:
            gene_annot[r["gene_id"]] = r["combined_name"]

print(f"Read in a total of {len(gene_annot):,} annotations")

fpi = "clustered.genes.fasta.gz"
fpo = "centroids.annot.csv.gz"

# Open both file paths, input and output
with gzip.open(fpi, "rt") as i, gzip.open(fpo, "wt") as o:

    # Write a header line
    o.write("gene_id,combined_name\\n")

    # Iterate over each input line
    for line in i:

        # If the line starts with '>'
        if line.startswith(">"):

            # Get the gene ID
            gene_id = line[1:].rstrip("\\n").split(" ")[0]

            # Write out the gene and any annotation that was found
            o.write(f"{gene_id},{gene_annot.get(gene_id, '')}\\n")
