#!/usr/bin/env python3
import logging
from pathlib import Path
import pandas as pd
import click

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [bin_summary.py] %(message)s'
)
logger = logging.getLogger('bin_summary.py')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Write to file
fileHandler = logging.FileHandler("bin_summary.log")
fileHandler.setFormatter(logFormatter)
fileHandler.setLevel(logging.INFO)
logger.addHandler(fileHandler)


class Counter:
    def __init__(self):
        pass

    def set_total(self, total: int):
        """
        Set the total number of items to process.
        """
        self.total = total
        self.current = 0
        self.reporting_interval = int(total / 100) if total > 100 else 1

    def increment(self, step: int = 1):
        """
        Increment the current count by a step.
        """
        self.current += step
        if self.current % self.reporting_interval == 0:
            logger.info(f"Processed {self.current:,} of {self.total:,} items.")


counter = Counter()


def main(read_alignments: pd.DataFrame, gene_bins: pd.DataFrame):
    """
    Main function to summarize the read alignments and gene bins.
    This function is called by the bin_metagenomes module.
    """
    logger.info("Calculating bin size - number of genes and total sequence length")
    bin_ngenes = gene_bins["bin"].value_counts().to_dict()
    bin_length_aa = gene_bins.groupby("bin")["length_aa"].sum().to_dict()

    # Get the total number of reads aligned to this pangenome per specimen
    tot_reads = read_alignments.groupby("specimen")["nreads"].sum()

    counter.set_total(read_alignments.reindex(["specimen", "bin"]).drop_duplicates().shape[0])

    logger.info("Summarizing read alignments by bin")
    (
        read_alignments.assign(bin=lambda d: d["id"].map(gene_bins["bin"].get))
        .groupby(["specimen", "bin"])
        .apply(lambda d: summarize_bin(d, bin_ngenes, bin_length_aa, tot_reads))
        .to_csv("bin_summary.csv.gz")
    )


def summarize_bin(df: pd.DataFrame, bin_ngenes: dict, bin_length_aa: dict, tot_reads_by_specimen: pd.Series) -> pd.Series:
    """
    Summarize the read alignments for a single bin.
    """
    # Number of genes with at least one read detected
    n_genes_detected = df["id"].nunique()

    # Bin ID
    bin_id = df["bin"].iloc[0]

    # Proportion of genes detected in this bin    
    prop_genes_detected = n_genes_detected / bin_ngenes[bin_id]

    # Coverage of genes in this bin, weighted for gene length
    bin_coverage = prop_genes_detected * (df["coverage"] * df["length_aa"]).sum() / df["length_aa"].sum()

    # Get the number of reads sequenced for this specimen
    tot_reads = tot_reads_by_specimen.loc[df["specimen"].values[0]]

    # Number and proportion of reads aligned
    n_reads_aligned = df["nreads"].sum()
    prop_reads_aligned = n_reads_aligned / tot_reads

    # Compute the RPKM (Reads Per Kilobase of transcript per Million mapped reads)
    rpkm = (
        (n_reads_aligned / (3 * bin_length_aa[bin_id] / 1000)) / (tot_reads / 1_000_000)
    )

    # Log the progress
    counter.increment()

    return pd.Series({
        "n_reads_aligned": n_reads_aligned,
        "prop_reads_aligned": prop_reads_aligned,
        "n_genes_detected": n_genes_detected,
        "prop_genes_detected": prop_genes_detected,
        "bin_size": bin_ngenes[bin_id],
        "bin_length_aa": bin_length_aa[bin_id],
        "rpkm": rpkm,
        "prop_genes_coverage": bin_coverage
    })


def read_csv(file_path: Path) -> pd.DataFrame:
    """
    Read a CSV file into a DataFrame.
    """
    logger.info(f"Reading in {file_path}")
    return pd.read_csv(file_path)


@click.command()
@click.option('--gene-bins', required=True, type=click.Path(exists=True), help='Path to gene_bins CSV file')
@click.option('--centroids-length', required=True, type=click.Path(exists=True), help='Path to centroids_length CSV file')
def cli(gene_bins, centroids_length):
    logger.info("Reading in all files within read_alignments/")
    read_alignments_df = pd.concat([
        read_csv(f)
        for f in Path("read_alignments").glob("*.csv.gz")
    ])

    logger.info(f"Reading in {gene_bins}")
    gene_bins_df = pd.read_csv(gene_bins)
    gene_bins_df.set_index("gene_id", inplace=True)

    centroids_length_df = pd.read_csv(centroids_length, index_col="header")["length"]

    gene_bins_df = gene_bins_df.assign(
        length_aa=lambda d: d.index.map(centroids_length_df.get)
    )

    read_alignments_df = read_alignments_df.assign(
        length_aa=lambda d: d["id"].map(centroids_length_df.get)
    )

    main(read_alignments_df, gene_bins_df)
    logger.info("Writing output to bin_summary.csv.gz")
    logger.info("Done.")


if __name__ == "__main__":
    cli()
