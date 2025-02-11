#!/usr/bin/env python3

import gzip
import os

input_fp = "${input_fasta}"
assert os.path.exists(input_fp)
output_fp = "${input_fasta.name}.filtered.fasta.gz"
print(f"Reading in {input_fp}, writing to {output_fp}")

min_gene_len = ${params.min_gene_length}
assert isinstance(min_gene_len, int), min_gene_len
print(f"Removing all genes shorter than {min_gene_len}aa")

def parse_fasta(handle):

    header = None
    seq = []
    counter = 0
    
    for line in handle:
        if line[0] == ">":
            if header is not None and len(seq) > 0:
                yield header, "".join(seq)
                counter += 1
            header = line[1:].split(" ")[0].rstrip("\\n")
            seq = []
        else:
            if len(line) > 1:
                seq.append(line.rstrip("\\n"))
                    
    if header is not None and len(seq) > 0:
        yield header, "".join(seq)
        counter += 1

    print(f"Read in a total of {counter:,} FASTA records")


def is_gzip_compressed(fp):
    with gzip.open(fp, 'rt') as handle:
        try:
            handle.read(1)
            return True
        except:
            return False


if is_gzip_compressed(input_fp):
    input_handle = gzip.open(input_fp, "rt")
else:
    input_handle = open(input_fp, "r")

counter = 0
with gzip.open(output_fp, 'wt') as output_handle:

    for header, seq in parse_fasta(input_handle):

        if len(seq) >= min_gene_len:

            output_handle.write(f">{header}\\n{seq}\\n")
            counter += 1

print(f"Wrote out {counter:,} sequences")

input_handle.close()
print("DONE")
