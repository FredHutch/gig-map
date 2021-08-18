#!/usr/bin/env python3
"""Select a set of marker genes from the aligned genes."""

import argparse
import logging
import os
import pandas as pd

def select_markers(
    alignments,
    min_coverage=50.,
    max_n=10
):

    # Calculate the coverage
    alignments = alignments.assign(
        coverage = 100 * alignments.length / alignments.slen
    )

    # Filter by coverage
    alignments = alignments.query(
        f"coverage >= {min_coverage}"
    )

    # Summarize each gene
    gene_df = alignments.groupby("sseqid").apply(summarize_gene)

    # Calculate the score for each gene
    gene_df = gene_df.assign(
        score = \
            gene_df.avg_iden * \
            (gene_df.n_genomes / gene_df.n_genomes.max()) * \
            (gene_df.gene_len / gene_df.gene_len.median()).clip(upper=1)
    ).sort_values(
        by=["score", "gene_len"],
        ascending=False
    )

    # Select the top N genes which are found in at least 3 genomes
    marker_genes = gene_df.query(
        "n_genomes >= 3"
    ).head(
        max_n
    ).index.values

    return marker_genes

def summarize_gene(df):

    # Get the number of genomes
    n_genomes = df.genome.unique().shape[0]

    # Get the average amino acid identity
    avg_iden = df.pident.mean() / 100.

    # Get the gene length
    gene_len = df.slen.values[0]
    
    return pd.Series(
        dict(
            n_genomes=n_genomes,
            avg_iden=avg_iden,
            gene_len=gene_len
        )
    )

if __name__ == "__main__":
    ##################
    # SET UP LOGGING #
    ##################

    # Set the level of the logger to INFO
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [select_markers.py] %(message)s'
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
        description="Select a set of marker genes from the aligned genes"
    )

    # Add the arguments
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Alignments of genes across genomes in CSV format'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Subset of alignments to just those selected marker genes'
    )
    parser.add_argument(
        '--max-n',
        type=int,
        required=True,
        help='Maximum number of marker genes to select'
    )
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=50.,
        help='Minimum alignment coverage threshold for selecting marker genes'
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
    alignments = pd.read_csv(args.input)

    # Filter down to the selected markers
    marker_genes = select_markers(
        alignments,
        min_coverage=args.min_coverage,
        max_n=args.max_n
    )

    # Write out the list of genes to a file
    with open(args.output, "w") as handle:

        handle.write("\n".join(marker_genes))
