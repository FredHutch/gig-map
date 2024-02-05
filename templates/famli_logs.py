#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import os

dat = defaultdict(dict)

for fp in os.listdir("."):
    if fp.endswith(".log"):
        sample = fp[:-len(".log")]
        with open(fp) as handle:
            for line in handle:
                if "[FAMLI parse] " not in line:
                    continue
                line = line.rstrip("\\n").split("[FAMLI parse] ", 1)[1]
                if line.startswith("Time elapsed: "):
                    dat[sample]["seconds"] = float(line[len("Time elapsed: "):].replace(",", ""))
                elif line.startswith("Results: assigned "):
                    line = line[len("Results: assigned "):]
                    assert line.endswith(" subjects")
                    line = line[:-len(" subjects")]
                    assert " queries to " in line
                    queries, subjects = line.split(" queries to ")
                    dat[sample]["queries"] = int(queries.replace(",", ""))
                    dat[sample]["subjects"] = int(subjects.replace(",", ""))
dat = pd.DataFrame(dat).T.reset_index().rename(
    columns=dict(
        index="sample"
    )
)
dat.to_csv("famli_logs.csv", index=None)