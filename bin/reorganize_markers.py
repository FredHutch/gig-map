#!/usr/bin/env python3

from collections import defaultdict
import gzip
import os

# Yet another FASTA parser
def read_fasta(handle):

    header = None
    seq = ""

    for line in handle:

        line = line.rstrip('\n')
        if line[0] == ">":

            if header is not None:
                yield header, seq

            header = line[1:].split(" ", 1)[0]
            seq = ""

        else:

            seq = f"{seq}{line}"

    yield header, seq


# The output object will be a dict of dicts
# The first-level key will be the marker of interest
# The second level key will be the genome of interest
# The final value will be the marker sequence for that gene in that genome
output = defaultdict(dict)


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
        output[marker_name][genome_name] = marker_sequence

# Iterate over the markers
for marker_name, marker_dict in output.items():

    # Open a file path for the output
    with gzip.open(
        os.path.join(output_folder, f"{marker_name}{output_file_ending}"),
        "wt"
    ) as handle:

        # Iterate over the sequences for each genome
        for genome_name, marker_sequence in marker_dict.items():

            # Write out the sequence
            handle.write(f">{genome_name}\n{marker_sequence}\n")
