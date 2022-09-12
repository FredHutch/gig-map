#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from statsmodels.stats.multitest import multipletests

# Get all of the long-form results
long_df = pd.concat([
    pd.read_csv(
        fp
    )
    for fp in Path(".").iterdir()
    if fp.name.endswith(".csv") and fp.name.startswith("corncob.results.")
])
long_df.to_csv("all.results.csv.gz", index=None)

# Make a set of wide-form results for each parameter
for parameter, df in long_df.groupby("parameter"):

    df = df.pivot_table(
        index="gene_id",
        columns="type",
        values="value"
    ).dropna(
    )

    if df.shape[0] == 0:
        continue

    fdr_bh = pd.Series(
        multipletests(df["p_value"].values, 0.1, method="fdr_bh")[1],
        index=df.index
    )
    
    df = df.assign(
        fdr_bh=fdr_bh
    )
    
    df.to_csv(
        f"{parameter.replace('(', '').replace(')', '')}.results.csv.gz"
    )
