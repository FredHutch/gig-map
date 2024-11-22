#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

(
    pd
    .concat([
        pd.read_csv(
            fp
        )
        for fp in Path(".").iterdir()
        if fp.name.endswith(".csv") or fp.name.endswith(".csv.gz")
    ])
    .to_csv("${params.output_csv}", index=None)
)
