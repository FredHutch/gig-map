#!/usr/bin/env python3

import gzip

fpi = "clustered.genes.fasta.gz"
fpo = "clustered.genes.csv.gz"

# Open both file paths, input and output
with gzip.open(fpi, "rt") as i, gzip.open(fpo, "wt") as o:

    # Write a header line
    o.write("gene_id,combined_name\\n")

    # Write to the output
    o.write(
        # A newline delimited list
        "\\n".join(
            [
                # Where each value is formatted from a line of the input
                # The formatting will 
                    # - remove the leading '>' character
                    # - replace the first ' ' with ','
                    # - and remove the trailing newline, if any
                line[1:].rstrip("\\n").replace(",", " ").replace(" ", ",", 1)
                # Parsed from each line of the input
                for line in i
                # For that subset of lines which start with '>' and which contain a space
                if line[0] == '>' and ' ' in line
            ]
        )
    )
