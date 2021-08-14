#!/usr/bin/env python3
"""Subset a table of alignments to just the specified genes."""

import argparse
import logging
import os
import pandas as pd


if __name__ == "__main__":
    ##################
    # SET UP LOGGING #
    ##################

    # Set the level of the logger to INFO
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [subset_alignments_by_genes.py] %(message)s'
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
        description="Subset a table of alignments to just the specified genes"
    )

    # Add the arguments
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Alignments of genes across genomes in TSV format'
    )
    # Add the arguments
    parser.add_argument(
        '--query-list',
        type=str,
        required=True,
        help='Text file with one gene name per line'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Subset of alignments to just those selected genes'
    )
    parser.add_argument(
        '--aln-fmt',
        type=str,
        default="qseqid sseqid pident length qstart qend qlen sstart send slen",
        help='Headers for the columns in the inputs, space delimited'
    )

    # Parse the arguments
    args = parser.parse_args()

    # Do not overwrite the output
    msg = f"Output already exists, stopping ({args.output}"
    assert os.path.exists(args.output) is False, msg

    # Read in the alignments
    assert os.path.exists(args.input)
    logger.info(f"Reading from {args.input}")
    alignments = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        names=args.aln_fmt.split(" ")
    )

    # Read in the query list
    assert os.path.exists(args.query_list)
    logger.info(f"Reading from {args.query_list}")
    query_list = open(args.query_list, 'r').readlines()

    query_list = [
        query_name.rstrip("\n")
        for query_name in query_list
        if len(query_name) > 1
    ]

    logger.info("Query list:")
    for i in query_list:
        logger.info(i)

    # Subset the alignments
    logger.info(f"Before filtering: {alignments.shape[0]:,} alignments")
    alignments = alignments.loc[
        alignments.sseqid.isin(set(query_list))
    ]
    logger.info(f"After filtering: {alignments.shape[0]:,} alignments")

    # Write out the filtered alignments to a file
    logger.info(f"Writing out to {args.output}")
    alignments.to_csv(
        args.output,
        sep="\t",
        header=None,
        index=None
    )
    logger.info("Done")