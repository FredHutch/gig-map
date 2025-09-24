#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("bin_summary.csv.gz")

for kw in ["rpkm", "n_reads_aligned", "prop_reads_aligned"]:
    assert kw in df.columns.values, f"Expected to find {kw} in the columns"
    print(f"Making wide table for {kw}")
    wide = df.pivot_table(
        index="specimen",
        columns="bin",
        values=kw
    ).fillna(
        0
    )
    print(f"Writing out table for {kw}")
    wide.to_csv(f"{kw}.csv.gz")

print("Done")
