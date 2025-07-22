#!/usr/bin/env python3
import gzip
import pandas as pd
from pathlib import Path
from collections import defaultdict


def compute_centroid_lengths(centroids_faa: Path) -> pd.DataFrame:
    lengths = defaultdict(int)
    header = None

    with gzip.open(centroids_faa, "rt") as f:
        for line in f:
            if line.startswith(">"):
                header = (
                    line[1:]
                    .strip()
                    .split(" ")[0]
                    .split("\t")[0]  # Use the first part of the header
                )
            else:
                lengths[header] += len(line.strip())
    return pd.DataFrame(lengths.items(), columns=["header", "length"])


def main(centroids_faa: Path, output_file: Path):
    """
    Main function to compute the lengths of centroids from a FASTA file.
    The results are saved to a specified output file.
    """
    lengths_df = compute_centroid_lengths(centroids_faa)
    lengths_df.to_csv(output_file, index=False, compression='gzip')
    print(f"Centroid lengths saved to {output_file}")


if __name__ == "__main__":
    main("${centroids_faa}", "centroids_length.csv.gz")
    # The placeholders \${centroids_faa} and \${centroids_length} will be replaced
    # by the actual file paths when this script is executed in the workflow.
