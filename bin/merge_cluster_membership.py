#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import gzip


def read_membership(fp: Path) -> pd.DataFrame:
    dat = []
    cluster_id = None
    with gzip.open(fp, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                cluster_id = line[1:].strip()
            else:
                member = line.strip().split(">")[1].split(" ")[0].strip(".")
                is_rep = line.strip().endswith("*")
                dat.append((cluster_id, member, is_rep))
    df = pd.DataFrame(dat, columns=["cluster", "member", "is_rep"])
    # Name each cluster for its representative
    rep_dict = df.query("is_rep").set_index("cluster")["member"].to_dict()
    df = df.assign(cluster=df["cluster"].map(rep_dict))
    return df.drop(columns=["is_rep"])


# Read in each of the cluster membership files from the shards
shards = pd.concat([
    read_membership(file)
    for file in Path(".").rglob("scatter.membership.*.tsv.gz")
])

# Read in the cluster membership file from the final gather across shards
gather = read_membership("gather.membership.tsv.gz")

# Replace the cluster ID from the shards with the cluster ID from the gather
shards = shards.assign(
    cluster=shards["cluster"].map(gather.set_index("member")["cluster"])
)

# Write out the merged cluster membership
shards.to_csv("centroids.membership.csv.gz", index=False)
