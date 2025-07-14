#!/usr/bin/env python3

from collections import defaultdict
import gzip
import sys
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Read in all FASTA files with the glob provided by sys.argv[0]
# Write out deduplicated names to a file with the prefix given by sys.argv[1]
# Drop any genes which are exactly the same
# If we encounter the same name twice, we will append a number to the end of the name
input_glob = sys.argv[1]
output_prefix = sys.argv[2]
output_suffix = "fasta.gz"
# Split up the output files with the maximum size cluster_shard_size
logger.info(f"Maximum shard size: {sys.argv[3]}")
cluster_shard_size = int(sys.argv[3])

logger.info(f"Reading in all FASTA files matching {input_glob}")
assert input_glob.endswith(".gz")
logger.info(f"Writing out deduplicated names to {output_prefix}.*.{output_suffix}")


def read_inputs():
    for fasta_fp in Path(".").rglob(input_glob):
        logger.info(f"Reading in {fasta_fp}")
        with gzip.open(fasta_fp, "rt") as handle:
            name = None
            seq = []
            for line in handle:
                if line.startswith(">"):
                    if name is not None and len(seq) > 0:
                        yield (name, ''.join(seq))
                    seq = []
                    name = line[1:].split(" ")[0].split("\t")[0].rstrip("\n")
                else:
                    seq.append(line.rstrip("\n"))
            if name is not None and len(seq) > 0:
                yield (name, ''.join(seq))


def drop_duplicates(input_generator):
    seen = set()
    for i in input_generator:
        if i not in seen:
            yield i
        seen.add(i)


all_names = defaultdict(int)
shard_ix = 1
shard_counter = 0

fpo = f"{output_prefix}.{shard_ix}.{output_suffix}"
logger.info(f"Writing out to {fpo}")
ofh = gzip.open(fpo, "wt")

for name, seq in drop_duplicates(read_inputs()):
    shard_counter += 1
    if shard_counter > cluster_shard_size:
        ofh.close()
        shard_ix += 1
        shard_counter = 0
        fpo = f"{output_prefix}.{shard_ix}.{output_suffix}"
        logger.info(f"Writing out to {fpo}")
        ofh = gzip.open(fpo, "wt")

    all_names[name] += 1
    if all_names[name] > 1:
        # If we have seen this name before, append a number to the end of the name
        ofh.write(f">{name}_{all_names[name]}\n{seq}\n")
    else:
        ofh.write(f">{name}\n{seq}\n")
ofh.close()

logger.info("Done")
