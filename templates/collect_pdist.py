#!/usr/bin/env python3

from collections import defaultdict
import gzip
import json
import os

# Both the input and output are formatted as
# {
#     min_name: {
#         max_name: pident
#     }
# }


def read_pdist(fp):
    with gzip.open(fp, 'rt') as handle:
        dat = json.load(handle)
    if len(dat) > 0:
        for min_name, dat in dat.items():
            for max_name, pident in dat.items():
                yield min_name, max_name, pident


output = defaultdict(dict)

# Iterate over each of the input files
for input_fp in os.listdir("."):
    if not input_fp.startswith("gene_pdist."):
        continue
    for min_name, max_name, pident in read_pdist(input_fp):
        output[min_name][max_name] = pident

# Write out the data
print("Writing out data")
with gzip.open("gene_pdist.json.gz", "wt") as handle:
    json.dump(output, handle)
print("Done")
