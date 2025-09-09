#!/usr/bin/env python3

import os
import pandas as pd
import gzip
import logging

logger = logging.getLogger('gig-map')
logger.setLevel(logging.INFO)
# Write to STDOUT
consoleHandler = logging.StreamHandler()
logger.addHandler(consoleHandler)


def read_gene_bins() -> pd.DataFrame:

    # Read in the gene bins file in CSV format
    bins_fp = "${gene_bins}"
    logger.info(f"Reading gene bins from {bins_fp}")
    bins = pd.read_csv(bins_fp)
    for cname in ['bin', 'gene_id']:
        if cname not in bins.columns:
            raise ValueError(f"Column '{cname}' not found in gene bins file")
    logger.info(f"Read in {len(bins):,} gene bins")
    return bins


def merge_bins():

    bins = read_gene_bins()

    # For each bin, make a single file with all of the distance matrix information,
    # and also merge the MSAs
    for bin_id, bin_df in bins.groupby('bin'):

        gene_list = sorted(bin_df['gene_id'].unique())

        logger.info(f"Merging data for bin {bin_id} with {len(gene_list):,} genes")

        # Merge the distance matrices
        merge_dms(bin_id, gene_list)
        merge_msas(bin_id, gene_list)
        logger.info(f"Finished merging data for bin {bin_id}")


def gene_distmat_fp(gene_id):
    return f"distmats/{gene_id}.distmat.csv.gz"


def flatten_dm(dm, gene_id):
    return (
        dm
        .reset_index()
        .melt(id_vars="index", var_name='genome2', value_name=gene_id)
        .set_index(['index', 'genome2'])
        .rename(index=dict(index='genome1'))
        [gene_id]
    )


def merge_dms(bin_id, gene_ids):
    # Read in each of the distance matrices
    dms = {
        gene_id: pd.read_csv(gene_distmat_fp(gene_id), index_col=0, compression="gzip")
        for gene_id in gene_ids
        if os.path.exists(gene_distmat_fp(gene_id))
    }
    if len(dms) == 0:
        return

    logger.info(f"Merging {len(dms):,} distance matrices for bin {bin_id}")

    # Flatten the long format distance matrices into a single DataFrame
    dms = pd.DataFrame({
        gene_id: flatten_dm(gene_df, gene_id)
        for gene_id, gene_df in dms.items()
    })

    dms.to_csv(f"{bin_id}.distmat.csv.gz", compression="gzip")


def merge_msas(bin_id, gene_ids):
    # Read each of the aligned sequences
    msas = {
        gene_id: read_fasta(f"msas/{gene_id}.msa.gz")
        for gene_id in gene_ids
        if os.path.exists(f"msas/{gene_id}.msa.gz")
    }
    if len(msas) == 0:
        return

    logger.info(f"Merging {len(msas):,} MSAs for bin {bin_id}")

    # For each alignment, make sure that the sequences are all the same length
    # Also keep track of the length of each gene
    gene_lengths = dict()
    for gene_id, msa in msas.items():
        lengths = set([len(seq) for seq in msa.values()])
        if len(lengths) != 1:
            raise ValueError(f"Sequences in MSA for gene {gene_id} are not all the same length")
        gene_lengths[gene_id] = lengths.pop()
        logger.info(f"Gene {gene_id} has {len(msa):,} sequences of length {gene_lengths[gene_id]:,}")

    # Get all of the unique genome IDs
    genomes = set([genome for gene in msas.values() for genome in gene.keys()])
    logger.info(f"Total of {len(genomes):,} unique genomes across all genes in bin {bin_id}")

    # Open up the output file which will have the combined MSA
    with gzip.open(f"{bin_id}.msa.gz", 'wt') as out_handle:

        # For each genome
        for genome_id in sorted(genomes):

            # Write the header
            out_handle.write(f">{genome_id}\\n")

            # For each gene, write the sequence if it exists, otherwise write Ns
            # (also replace any gaps with Ns)
            for gene_id, gene_msa in msas.items():

                if genome_id in gene_msa:
                    out_handle.write(gene_msa[genome_id].replace("-", "N"))
                else:
                    out_handle.write('N' * gene_lengths[gene_id])

            out_handle.write("\\n")

            logger.info(f"Wrote combined MSA for genome {genome_id}")


def read_fasta(fp):
    """
    Read a FASTA file and return a dictionary of sequences.

    Parameters:
    fp (str): Path to the input FASTA file.

    Returns:
    dict: A dictionary where keys are sequence headers and values are sequences.
    """

    header = None
    seq = []
    seqs = dict()

    with gzip.open(fp, 'rt') as handle:
        for line in handle:

            line = line.rstrip('\\n')

            # If the line is empty
            if len(line) == 0:

                # Skip it
                continue

            if line[0] == ">":

                if header is not None:
                    seqs[header] = ''.join(seq)

                header = line[1:].split(" ", 1)[0]
                seq = []

            else:

                seq.append(line)

    seqs[header] = ''.join(seq)

    return seqs


if __name__ == "__main__":
    merge_bins()
