#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import os

pattern = [
    {
        "prefix": "#CPU threads: ",
        "suffix": "",
        "key": "cpus",
        "type": int
    },
    {
        "prefix": "Total time = ",
        "suffix": "s",
        "key": "seconds",
        "type": float
    },
    {
        "prefix": "",
        "suffix": " queries aligned.",
        "key": "queries",
        "type": int
    }
]

dat = defaultdict(dict)

for fp in os.listdir("."):
    if fp.endswith(".log"):
        sample, index = fp[:-len(".log")].rsplit(".", 1)
        with open(fp) as handle:
            for line in handle:
                line = line.rstrip("\\n")
                for p in pattern:
                    if line.startswith(p["prefix"]) and line.endswith(p["suffix"]):
                        if len(p["prefix"]) > 0:
                            line = line[len(p["prefix"]):]
                        if len(p["suffix"]) > 0:
                            line = line[:-len(p["suffix"])]
                        val = p["type"](line)
                        dat[(sample, index)][p["key"]] = val
                        break
                    
dat = pd.DataFrame(dat).T.reset_index().rename(
    columns=dict(
        level_0="sample",
        level_1="index",
    )
)
dat.to_csv("alignment_logs.csv", index=None)