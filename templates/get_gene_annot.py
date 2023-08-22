#!/usr/bin/env python3

import gzip

fpi = "${fasta}"
fpo = "${fasta}.annot.csv.gz"


def safe_open(fp):
    if fp.endswith(".gz"):
        return gzip.open(fp, "rt")
    else:
        return open(fp, "r")


# Open both file paths, input and output
with safe_open(fpi) as i, gzip.open(fpo, "wt") as o:

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
                # For that subset of lines which start with '>'
                if line[0] == '>'
            ]
        )
    )
