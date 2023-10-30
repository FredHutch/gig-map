#!/usr/bin/env python3
"""Select a set of marker genes from the aligned genes."""

import argparse
import logging
import os
import gzip
import pandas as pd


def select_markers(
    alignments: pd.DataFrame,
    min_coverage=50.,
    max_n=10
):
    logger = logging.getLogger('gig-map')

    # Calculate the coverage
    alignments = alignments.assign(
        coverage=100 * alignments.length / alignments.slen
    )

    # Filter by coverage
    alignments = alignments.query(
        f"coverage >= {min_coverage}"
    )
    logger.info(f"Rows after filering by coverage: {alignments.shape[0]:,}")

    # Summarize each gene
    gene_df = alignments.groupby("sseqid").apply(summarize_gene)
    logger.info(f"Gene Summary:\n{gene_df.head().to_csv()}")

    # Calculate the score for each gene
    gene_df = gene_df.assign(
        score=calc_score(gene_df)
    ).sort_values(
        by=["score", "gene_len"],
        ascending=False
    )
    logger.info(f"Scored & Sorted:\n{gene_df.head().to_csv()}")

    # Select the top N genes which are found in at least 3 genomes
    marker_genes = gene_df.query(
        "n_genomes >= 3"
    ).head(
        max_n
    ).index.values

    return marker_genes


def calc_score(gene_df: pd.Series) -> float:
    logger = logging.getLogger('gig-map')
    n_max = gene_df.n_genomes.max()
    logger.info(f"Maximum number of genomes aligned: {n_max}")
    gl_m = gene_df.gene_len.median()
    logger.info(f"Median gene length: {gl_m}")
    return (
        gene_df
        .apply(
            lambda r: r.avg_iden * (r.n_genomes / n_max) * (r.gene_len / gl_m),
            axis=1
        )
        .clip(upper=1)
    )


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


def subset_fasta(gene_list, fasta_in, fasta_out):
    """Filter a FASTA file to just those sequences in the provided list"""

    gene_set = set(gene_list)

    # Open the input and output file handles
    with gzip.open(fasta_out, "wt") as handle_o, gzip.open(fasta_in, "rt") as handle_i:

        for header, seq in fasta_parser(handle_i):

            if header in gene_set:

                handle_o.write(f">{header}\n{seq}\n")


def fasta_parser(handle):
    """Parse a file handle in FASTA format."""

    header = None
    seq = []

    for line in handle:
        if line[0] == ">":
            if header is not None and len(seq) > 0:
                yield header, "".join(seq)
            header = line[1:].split(" ")[0].rstrip("\n")
            seq = []
        else:
            if len(line) > 1:
                seq.append(line.rstrip("\n"))

    if header is not None and len(seq) > 0:
        yield header, "".join(seq)


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
        '--all-genes',
        type=str,
        required=True,
        help='File containing all gene sequences in FASTA format (gzip-compressed)'
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
    logger.info(f"Read in {alignments.shape[0]:,} rows")

    # Filter down to the selected markers
    marker_genes = select_markers(
        alignments,
        min_coverage=args.min_coverage,
        max_n=args.max_n
    )

    # Write out the sequences of the selected markers
    subset_fasta(
        marker_genes,
        args.all_genes,
        args.output
    )
