process assign_metadata {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
    path input_metadata_csv

    output:
    path "metadata.csv"

    """#!/usr/bin/env python3

import pandas as pd
import numpy as np

print("Reading in ${input_metadata_csv}")
df = pd.read_csv("${input_metadata_csv}")
print(f"Read in {df.shape[0]:,} rows and {df.shape[1]:,} columns")

# If the assign_metadata parameter was set, use that
assign_metadata_str = "${params.assign_metadata}"
if len(assign_metadata_str) > 0 and assign_metadata_str != "false":
    print("Applying assignment: ${params.assign_metadata}")
    df = df.assign(${params.assign_metadata})


df.to_csv("metadata.csv", index=None)
print(df.shape)
print("Wrote out metadata.csv")
    """

}