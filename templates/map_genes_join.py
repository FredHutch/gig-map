#!/usr/bin/env python3

import os
import pandas as pd

combined = []
for fp in os.listdir("."):
    if fp.startswith("gene_mapping.shard"):
        print(f"Reading {fp}")
        df = pd.read_csv(fp)
        combined.append(df)
print(f"Joining {len(combined):,} tables")
combined = pd.concat(combined)

print("Writing out")

combined.to_csv("gene_mapping.csv.gz", sep=",", index=None)
