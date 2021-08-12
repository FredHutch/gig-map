#!/usr/bin/env python3

import gzip
import logging
import os
import pandas as pd
import sys

# Input arguments are <alignments> and <genome_fasta>
alignments_fp = sys.argv[1]
genome_fp = sys.argv[2]
header_string = sys.argv[3]
min_coverage = float(sys.argv[4])
assert os.path.exists(alignments_fp)
assert os.path.exists(genome_fp)

# Yet another FASTA parser
def read_fasta(handle):

    header = None
    seq = ""

    genome = dict()

    for line in handle:

        line = line.rstrip('\n')
        if line[0] == ">":

            if header is not None:
                genome[header] = seq

            header = line[1:].split(" ", 1)[0]
            seq = ""

        else:

            seq = f"{seq}{line}"

    genome[header] = seq

    return genome


# Function to extract the aligned portion of the genome
def get_aligned_region(aln_r, genome):

    # Make sure that we have the right contig in the FASTA
    assert aln_r.qseqid in genome.keys(), f"Cannot find {aln_r.qseqid} in {genome_fp}"

    # Make sure that the nucleotide sequence is long enough
    assert len(genome[aln_r.qseqid]) >= aln_r.qstart
    assert len(genome[aln_r.qseqid]) >= aln_r.qend

    # If the alignment is in the forward direction
    if aln_r.qstart < aln_r.qend:

        # Extract the region
        region = genome[aln_r.qseqid][
            int(aln_r.qstart - 1):int(aln_r.qend)
        ]

    # If the alignment is in the reverse direction
    else:

        # Extract the region
        region = genome[aln_r.qseqid][
            int(aln_r.qend - 1):int(aln_r.qstart)
        ]

        # Reverse complement it
        region = "".join(
            map(
                lambda nuc: dict(A="T", T="A", C="G", G="C", N="N").get(nuc, "N"),
                list(region.upper())[::-1]
            )
        )

    return region


# Read in the alignments
logging.info(f"Reading in {alignments_fp}")
aln = pd.read_csv(
    alignments_fp,
    sep="\t",
    header=None,
    names=header_string.split(" ")
)
logging.info(f"Read in {aln.shape[0]:,} alignments")

# Read in the genome

# If the file is gzipped
if genome_fp.endswith(".gz"):

    # Read it in using gzip
    handle = gzip.open(genome_fp, "rt")

# Otherwise
else:

    # Read it in uncompressed
    handle = open(genome_fp, "r")

genome = read_fasta(handle)
handle.close()
logging.info(f"Read in {len(genome):,} genome records")

# If there are alignments
if aln.shape[0] > 0:

    # Filter by the minimum coverage threshold
    aln = aln.assign(
        coverage = 100. * aln.length / aln.slen
    ).query(
        f"coverage >= {min_coverage}"
    )
    
    logging.info(f"{aln.shape[0]:,} alignments pass the coverage threshold of {min_coverage}")

    # If there are alignments which pass that threshold
    if aln.shape[0] > 0:

        # Assign a score based on alignment length and sequence identity
        aln = aln.assign(
            score = aln.pident * aln.length
        )

        # Make a dict with the aligned sequence for each marker
        marker_sequences = dict()

        # Iterate over each marker
        for marker_name, marker_aln in aln.groupby("sseqid"):

            # Get the aligned sequence for the top-scoring alignment
            marker_sequences[
                marker_name
            ] = get_aligned_region(
                marker_aln.sort_values(by="score", ascending=False).head(1).iloc[0],
                genome
            )

        # Write out to a file
        with gzip.open(f"{genome_fp}.markers.fasta.gz", "wt") as handle:
            for marker_name, marker_sequence in marker_sequences.items():
                handle.write(f">{marker_name}\n{marker_sequence}\n")