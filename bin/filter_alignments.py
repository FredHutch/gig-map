#!/usr/bin/env python3
"""Filter a set of alignments to remove overlapping alignment regions."""

import argparse
from bitarray import bitarray
import logging
import os
import pandas as pd
from time import sleep

def filter_alignments(alignments: pd.DataFrame, max_overlap: int):

    # Get the logger
    logger = logging.getLogger('gig-map')

    # While `max_overlap` is a percentage, we will use a proportion internally
    assert max_overlap >= 0
    assert max_overlap <= 100
    max_overlap = max_overlap / 100.
    assert isinstance(max_overlap, float)
    assert max_overlap > 0.
    assert max_overlap <= 1.

    # Assign a score to each alignment which is based on the coverage,
    # percent identity, and gene length

    # The gene lengths will be weighted as min(L / median(L), 1)

    # Find the median length of genes which are aligned to this genome
    logger.info("Calculating median gene length")
    median_gene_length = alignments.reindex(
        columns=["sseqid", "slen"]
    ).drop_duplicates(
    )[
        "slen"
    ].median()
    logger.info(f"Median gene length is {median_gene_length:,}")

    # Calculate the proportion of the gene which is aligned
    logger.info("Calculating alignment coverage")
    alignments = alignments.assign(
        coverage=(alignments.length / alignments.slen).clip(upper=1)
    )

    # Calculate the score per alignment
    logger.info("Calculating alignment scores")
    alignments = alignments.assign(
        aln_score=alignments.coverage * alignments.pident / 100. * (alignments.slen / median_gene_length).clip(upper=1)
    )

    # Sort by score
    logger.info("Sorting alignments by score")
    alignments = alignments.sort_values(
        by="aln_score",
        ascending=False
    ).reset_index(
        drop=True
    )

    # Keep a list of the index positions to keep
    alignments_to_keep = []

    # Make a dict with an empty bitarray for each contig
    logger.info("Initializing coverage arrays per contig")
    contig_coverage = {
        contig_name: empty_array(contig_len)
        for contig_name, contig_len in alignments.reindex(
            columns=["qseqid", "qlen"]
        ).drop_duplicates(
        ).set_index(
            "qseqid"
        )[
            "qlen"
        ].items()
    }

    # Iterate over each alignment
    for i, r in alignments.iterrows():

        # print(i)
        # print(r)

        # Calculate the overlap of this region with previously aligned
        # genes for this particular contig
        overlap_prop = calc_overlap(r, contig_coverage[r.qseqid])

        # If the overlap is under the threshold
        if overlap_prop < max_overlap:

            # Add it to the list
            alignments_to_keep.append(i)

            # Mark the region as being covered
            add_alignment(r, contig_coverage[r.qseqid])

    logger.info(f"Keeping {len(alignments_to_keep):,} / {alignments.shape[0]:,} alignments")

    # Filter the table and sort it
    alignments = alignments.reindex(
        index=alignments_to_keep
    ).sort_values(
        by=["qseqid", "qstart", "qend"]
    )

    # Delete the added columns
    alignments = alignments.drop(
        columns=["aln_score", "coverage"]
    )

    return alignments


def find_aln_start_stop(r):
    """
    Find the 5' and 3' positions for the alignment
    """

    # Make sure to adjust for zero-indexing
    if r.qstart < r.qend:
        aln_left = r.qstart - 1
        aln_right = r.qend
    else:
        aln_left = r.qend - 1
        aln_right = r.qstart

    return aln_left, aln_right


def add_alignment(r, prev_map):
    """Mark the region as being covered."""

    # Find the 5' and 3' positions for the alignment
    aln_left, aln_right = find_aln_start_stop(r)

    # Set the bits to 1
    prev_map[aln_left:aln_right] = 1


def calc_overlap(r, prev_map):
    """Calculate the proportion of this alignment which is already covered."""

    # Find the 5' and 3' positions for the alignment
    aln_left, aln_right = find_aln_start_stop(r)

    # Get the number of bases which overlap
    n_overlap = prev_map.count(1, aln_left, aln_right)

    # Divide by the length of the alignment
    prop_overlap = n_overlap / (aln_right - aln_left)

    return prop_overlap


def empty_array(contig_len: int):
    """Make an empty bitarray of a given length."""
    a = bitarray(contig_len)
    a.setall(0)
    return a


if __name__ == "__main__":

    ##################
    # SET UP LOGGING #
    ##################

    # Set the level of the logger to INFO
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [filter_alignments.py] %(message)s'
    )
    logger = logging.getLogger('gig-map')
    logger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    ###################
    # PARSE ARGUMENTS #
    ###################

    # Create the parser
    parser = argparse.ArgumentParser(
        description="Filter a set of alignments to remove overlapping alignment regions"
    )

    # Add the arguments
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Alignments of genes across genomes in TSV format'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Filtered alignments of genes across genomes in TSV format'
    )
    parser.add_argument(
        '--max-overlap',
        type=int,
        required=True,
        help='Remove any alignments which overlap a higher scoring alignment by at least this percentage (0 - 100)'
    )
    parser.add_argument(
        '--aln-fmt',
        type=str,
        default="qseqid sseqid pident length qstart qend qlen sstart send slen",
        help='Headers for the columns in the inputs, space delimited'
    )

    # Parse the arguments
    args = parser.parse_args()

    ##########
    # INPUTS #
    ##########

    # Make sure that the output file doesn't already exist
    msg = f"{args.output} already exists -- stopping"
    assert os.path.exists(args.output) is False, msg

    # Make sure that the --max-overlap falls in the range 0-100
    msg = "--max-overlap must be an integer in the range 0 - 100"
    assert isinstance(args.max_overlap, int), msg
    assert args.max_overlap >= 0, msg
    assert args.max_overlap <= 1000, msg

    # Read in the alignments
    logger.info(f"Reading from {args.input}")
    assert os.path.exists(args.input)
    alignments = pd.read_csv(
        args.input,
        header=None,
        names=args.aln_fmt.split(" "),
        sep="\t"
    )

    ##########
    # FILTER #
    ##########

    # Filter the alignments
    logger.info(f"Filtering alignments with max overlap {args.max_overlap}")
    filtered_alignments = filter_alignments(
        alignments,
        args.max_overlap
    )

    ##########
    # OUTPUT #
    ##########

    # Write out to a file
    logger.info(f"Writing out to {args.output}")
    filtered_alignments.to_csv(
        args.output,
        index=None,
        header=None,
        sep="\t"
    )
