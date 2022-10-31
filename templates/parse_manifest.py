#!/usr/bin/env python3
"""Make sure that the manifest provided by the user is valid."""

import pandas as pd

print("Reading in the CSV file specified by --manifest (${params.manifest})")
manifest = pd.read_csv("raw.manifest.csv")
print(f"Read in {manifest.shape[0]:,} lines from ${params.manifest}")
assert "specimen" in manifest.columns.values, "ERROR: manifest must contain 'specimen' column"

print("Reading in the CSV file specified by --gene_abund (${params.gene_abund})")
gene_abund = pd.read_csv("read_alignments.csv.gz", usecols=["specimen", "id"])
print(f"Read in {gene_abund.shape[0]:,} lines from ${params.gene_abund}")
assert "specimen" in gene_abund.columns.values, "ERROR: gene abundance CSV must contain 'specimen' column"
assert "id" in gene_abund.columns.values, "ERROR: gene abundance CSV must contain 'id' column"

# Get the set of specimens in the manifest and abundance tables
manifest_specimens = set(manifest["specimen"].values)
print(f"Specimens in the manifest: {', '.join(list(manifest_specimens))}")
gene_abund_specimens = set(gene_abund["specimen"].values)
print(f"Specimens in the gene abundance CSV: {', '.join(list(gene_abund_specimens))}")

# All of the specimens in the manifest must be in the gene abundances
missing_specimens = manifest_specimens - gene_abund_specimens
msg = f"ERROR: All specimens in the manifest must be in the gene abundance CSV ({', '.join(list(missing_specimens))})"
assert len(missing_specimens) == 0, msg

# Get the columns listed in the formula
formula = "${params.formula}"
formula_cols = formula.split(" + ")
missing_cols = set(formula_cols) - set(manifest.columns.values)
msg = f"ERROR: Columns specified in the formula are not present in the manifest: {', '.join(list(missing_cols))}"
assert len(missing_cols) == 0, msg

# Subset to the formula columns and drop any missing values
manifest = manifest.set_index(
    "specimen"
).reindex(
    columns=formula_cols
).dropna(
).reset_index(
).drop_duplicates()
print(f"Rows with valid data for {', '.join(formula_cols)}: {manifest.shape[0]:,}")

# Write out the manifest
manifest.to_csv(
    "parsed.manifest.csv"
)