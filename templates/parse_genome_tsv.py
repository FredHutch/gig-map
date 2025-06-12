#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("input.tsv", sep="\\t")
msg = "Did not find expected column: Assembly Accession"
assert "Assembly Accession" in df.columns.values, msg

# Make a list of accessions to download
acc_list = []


# Iterate over each row in the table
for _, r in df.iterrows():

    acc = r["Assembly Accession"]

    if (
        (not isinstance(acc, str)) or
        (len(acc) == 0)
    ):
        print("No valid assembly accession found for record:")
        print(r.to_json())

    if acc in acc_list:
        raise ValueError(f"Duplicate assembly accession found: {acc}")

    acc_list.append(acc)


# Write the list to a file
with open("acc_list.txt", "w") as handle:

    # Each accession on its own line
    handle.write("\\n".join(acc_list))
