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

    print(f"Reading in {input_fp}")

    # Read in the sequence for each marker
    # The header for each record is:
    #     >{marker_name}::{genome_name}
    for marker_header, marker_sequence in read_fasta(
        os.path.join(
            input_folder, input_fp
        )
    ):

        # Make sure that the appropriate delimiter is found
        msg = f"Expected to find '::' in '{marker_header}'"
        assert "::" in marker_header, msg

        # Parse the marker name and genome name from the header
        marker_name, genome_name = marker_header.split("::", 1)

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
print(f"Read in {df.shape[0]:,} records in total")

# Iterate over the markers
for marker_name, marker_df in df.groupby("marker_name"):

    print(f"Reorganizing genes for {marker_name}")

    # Make sure that we only write out a single sequence for each genome
    # This should only be an edge case that comes up when the user provides
    # a marker gene which is also the one selected automatically
    written_genomes = set()

    # Open a file path for the output
    fpo = os.path.join(output_folder, f"{marker_name}{output_file_ending}")
    print(f"Writing out to {fpo}")
    with gzip.open(fpo, "wt") as handle:

        # Iterate over the sequences for each genome
        for _, r in marker_df.iterrows():

            # If we have already written out this genome
            if r.genome_name in written_genomes:

                # Skip it
                continue

            # If we have not yet written this genome
            else:

                # Write out the sequence
                handle.write(f">{r.genome_name}\n{r.marker_sequence}\n")

                # Note that we've written out this genome
                written_genomes.add(r.genome_name)
print("DONE")