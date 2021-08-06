#!/usr/bin/env python3

from collections import defaultdict
import gzip
import os
import pandas as pd

# Yet another FASTA parser
def read_fasta(fp):

    if fp.endswith(".gz"):
        lines = gzip.open(fp, "rt").readlines()
    else:
        lines = open(fp, "r").readlines()

    header = None
    seq = ""

    for line in lines:

        line = line.rstrip('\n')
        if line[0] == ">":

            if header is not None:
                yield header, seq

            header = line[1:].split(" ", 1)[0]
            seq = ""

        else:

            seq = f"{seq}{line}"

    yield header, seq


# Organize markers in a table
output = []


# Make variables with the input and output folders and file endings
input_folder = "fastas_by_genome/"
input_file_ending = ".markers.fasta.gz"
output_folder = "fastas_by_marker/"
output_file_ending = ".markers.fasta.gz"


# Walk through the files in fastas_by_genome/
for input_fp in os.listdir(input_folder):

    # Each file must have the suffix .markers.fasta.gz
    assert input_fp.endswith(input_file_ending)

    # Get the genome name from the file name
    genome_name = input_fp.replace(input_file_ending, "")

    # Read in the sequence for each marker
    for marker_name, marker_sequence in read_fasta(
        os.path.join(
            input_folder, input_fp
        )
    ):

        # Add it to the output
        output.append(
            dict(
                marker_name=marker_name,
                marker_sequence=marker_sequence,
                genome_name=genome_name
            )
        )

# Format as a DataFrame
df = pd.DataFrame(output)

# Iterate over the markers
for marker_name, marker_df in df.groupby("marker_name"):

    # Open a file path for the output
    with gzip.open(
        os.path.join(output_folder, f"{marker_name}{output_file_ending}"),
        "wt"
    ) as handle:

        # Iterate over the sequences for each genome
        for _, r in marker_df.iterrows():

            # Write out the sequence
            handle.write(f">{r.genome_name}\n{r.marker_sequence}\n")
