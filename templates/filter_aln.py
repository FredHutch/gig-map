#!/usr/bin/env python3

import gzip
import os

# File containing the alignments
aln = "${aln}"
assert os.path.exists(aln), f"Did not find expected file {aln}"

# Count the number of lines
n = 0
with gzip.open(aln, "rt") as handle:
    for line in handle:
        n += 1

print(f"Found {n:,} lines in {aln}")

# If this number is below the threshold
threshold = int("${params.min_alignments}")
if n < threshold:
    print(f"Falls below threshold ({threshold:,}), removing")

    # Delete the file
    os.remove(aln)

else:
    print(f"Satisfies threshold ({threshold:,})")

print("Done")
