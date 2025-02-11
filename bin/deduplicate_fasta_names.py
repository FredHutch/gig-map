#!/usr/bin/env python3

from collections import defaultdict
import gzip
import sys
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.INFO)

# Read in all FASTA files with the glob provided by sys.argv[0]
# Write out deduplicated names to the path given by sys.argv[1]
# If we encounter the same name twice, we will append a number to the end of the name
input_glob = sys.argv[1]
output_fp = sys.argv[2]
logger.info(f"Reading in all FASTA files matching {input_glob}")
logger.info(f"Writing out deduplicated names to {output_fp}")
assert output_fp.endswith(".gz")

all_names = defaultdict(int)

with gzip.open(output_fp, "wt") as ofh:
    for fasta_fp in Path(".").rglob(input_glob):
        logger.info(f"Reading in {fasta_fp}")
        assert fasta_fp.name.endswith(".gz")
        with gzip.open(fasta_fp, "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split(" ")[0].rstrip("\n")
                    all_names[name] += 1
                    if all_names[name] > 1:
                        logger.info(f"Encountered {name} {all_names[name]:,} times")
                        ofh.write(f">{name}_{all_names[name]}\n")
                    else:
                        ofh.write(line)
                else:
                    ofh.write(line)
logger.info("Done")
