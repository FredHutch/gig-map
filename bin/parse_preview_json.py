#!/usr/bin/env python3
"""Parse the preview.json produced by the NCBI datasets CLI"""

import json

dat = json.load(open('preview.json'))
print("Checking for number of records")
print(json.dumps(dat, indent=4))

for kw, val in dat["included_data_files"].items():
    if val["file_count"] > 0:
        with open("accession.has.data", "w") as handle:
            handle.write('TRUE')
        break
