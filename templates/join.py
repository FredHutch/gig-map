#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from statsmodels.stats.multitest import multipletests

# Get all of the long-form results
df = pd.concat([
    pd.read_csv(fp)
    for fp in Path(".").iterdir()
    if fp.name.endswith(".csv") and fp.name.startswith("regress.results.")
])

# Add an FDR-BH for everything together
df = df.assign(fdr_bh=multipletests(df["pvalue"], method="fdr_bh")[1])

df.to_csv("regress.results.csv.gz", index=None)