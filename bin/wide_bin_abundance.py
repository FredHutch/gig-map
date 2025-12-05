#!/usr/bin/env python3

import pandas as pd
import click


def filter_specimens(df: pd.DataFrame, val: int, kw: str) -> pd.DataFrame:
    # Get the sum of the kw column for each specimen
    specimen_vals = df.groupby("specimen")[kw].sum()

    # Get the specimens which pass the filter
    passing = specimen_vals.index.values[specimen_vals >= val]

    print(f"Specimens with >= {val} for {kw}: {len(passing):,} / {specimen_vals.shape[0]:,}")
    return df.loc[df["specimen"].isin(passing)]


@click.command()
@click.option('--min-n-reads', type=int, required=True, help='Minimum number of reads to include a sample')
@click.option('--min-n-genes', type=int, required=True, help='Minimum number of genes to include a sample')
def main(min_n_reads, min_n_genes):

    df = pd.read_csv("bin_summary.csv.gz")
    print(f"Minimum number of reads to include a sample: {min_n_reads:,}")
    print(f"Minimum number of genes to include a sample: {min_n_genes:,}")

    df = filter_specimens(df, min_n_reads, "n_reads_aligned")
    df = filter_specimens(df, min_n_genes, "n_genes_detected")

    for kw in ["rpkm", "n_reads_aligned", "prop_reads_aligned"]:
        assert kw in df.columns.values, f"Expected to find {kw} in the columns"
        print(f"Making wide table for {kw}")
        wide = df.pivot_table(
            index="specimen",
            columns="bin",
            values=kw
        ).fillna(
            0
        )
        print(f"Writing out table for {kw}")
        wide.to_csv(f"{kw}.csv.gz")

        # Use this to output the "fragments-per-million" table
        if kw == "prop_reads_aligned":
            print("Writing out the table for fragments_per_million")
            (wide * 1e9).to_csv("fragments_per_million.csv.gz")

    print("Done")

if __name__ == "__main__":
    main()
