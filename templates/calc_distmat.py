#!/usr/bin/env python3

import pandas as pd
import gzip


# Yet another FASTA parser
def read_fasta(fp):

    header = None
    seq = []

    genome = dict()

    with gzip.open(fp, 'rt') as handle:
        for line in handle:

            line = line.rstrip('\\n')

            # If the line is empty
            if len(line) == 0:

                # Skip it
                continue

            if line[0] == ">":

                if header is not None:
                    genome[header] = ''.join(seq)

                header = line[1:].split(" ", 1)[0]
                seq = []

            else:

                seq.append(line)

    genome[header] = ''.join(seq)

    return genome


def calc_distmat(msa_file, output_file):
    """
    Calculate a distance matrix from a multiple sequence alignment (MSA) file.
    Note - any missing positions will be ignored in the distance calculation.

    Parameters:
    msa_file (str): Path to the input MSA file in FASTA format.
    output_file (str): Path to the output CSV file to save the distance matrix.

    Returns:
    None
    """

    # Read the MSA file - which is in FASTA format
    sequences = read_fasta(msa_file)

    # Calculate the pairwise distance matrix
    dm = pd.DataFrame({
        seq_id1: {
            seq_id2: (
                calc_dm(seq1, seq2)
                if seq_id1 != seq_id2 else 0.0
            )
            for seq_id2, seq2 in sequences.items()            
        }
        for seq_id1, seq1 in sequences.items()
    }).sort_index().sort_index(axis=1)

    # Save the distance matrix to a CSV file
    dm.to_csv(output_file, index=True, header=True, float_format='%.6f')


def calc_dm(seq1, seq2):
    """
    Calculate the distance between two sequences as the proportion of differing positions,
    ignoring any positions with gaps or missing data.

    Parameters:
    seq1 (str): The first sequence.
    seq2 (str): The second sequence.

    Returns:
    float: The distance between the two sequences.
    """

    assert len(seq1) == len(seq2), "Sequences must be of equal length"

    differences = 0
    valid_positions = 0

    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        if a != 'N' and b != 'N':
            valid_positions += 1
            if a != b:
                differences += 1

    if valid_positions == 0:
        return 0.0

    return differences / valid_positions


if __name__ == "__main__":
    calc_distmat("${msa_fasta}", "${msa_fasta}".replace('.msa.gz', '.distmat.csv.gz'))
