#!/usr/bin/env python3

from typing import List
import pandas as pd
from pathlib import Path
import plotly.express as px
from scipy.cluster import hierarchy
import logging

logger = logging.getLogger(__name__)


def read_file(fp: Path, header: List[str]):
    logger.info(f"Reading data from {fp}")

    return pd.read_csv(
        fp,
        sep="\\t",
        header=None,
        names=header
    )


def run():

    # Get the alignment information
    aln = read_aln()

    # Plot the number of genomes each gene was aligned to
    gene_frequency_bars(aln)

    # Plot the number of genomes each gene was aligned to
    genome_frequency_bars(aln)

    # Plot the distribution of alignment coverage and identity metrics
    alignment_qc_metrics(aln)

    # Plot a heatmap of all aligned genes
    alignment_heatmap(aln)

    # Count up the number of genes per genome, given different thresholds
    gene_frequency_genome_size(aln)


def read_aln() -> pd.DataFrame:

    # Parse the expected column names
    header = "${params.aln_fmt} genome".split(" ")
    logger.info(f"Header: {','.join(header)}")

    # Read in all of the files and join the tables
    aln: pd.DataFrame = pd.concat([
        read_file(fp, header)
        for fp in Path("inputs").rglob("*.tsv.gz")
    ])

    # Add the coverage information
    aln = aln.assign(
        coverage=(100 * aln['length'] / aln['slen']).clip(upper=100.)
    )

    # Write out the concatenated alignments
    aln.to_csv("genomes.aln.csv.gz")

    ngenes = aln['sseqid'].unique().shape[0]
    ngenomes = aln['genome'].unique().shape[0]
    logger.info(f"Loaded {aln.shape[0]:,} alignments")
    logger.info(f"{ngenes:,} genes and {ngenomes:,} genomes")

    return aln


def gene_frequency_bars(aln: pd.DataFrame):

    # Count up the number of genomes each gene was aligned to
    vc = (
        aln
        .reindex(columns=["genome", "sseqid"])
        .drop_duplicates()
        ["sseqid"]
        .value_counts()
        .reset_index()
    )

    # Make a bargraph
    fig = px.bar(
        data_frame=vc,
        x="index",
        y="sseqid",
        labels=dict(sseqid="Number of Genomes", index="Gene"),
        title="Gene Detection Frequency"
    )
    fig.write_html("gene_frequency_bars.html", include_plotlyjs="cdn")


def genome_frequency_bars(aln: pd.DataFrame):

    # Count up the number of genes each genome contains
    vc = (
        aln
        .reindex(columns=["genome", "sseqid"])
        .drop_duplicates()
        ["genome"]
        .value_counts()
        .reset_index()
    )

    # Make a bargraph
    fig = px.bar(
        data_frame=vc,
        x="index",
        y="genome",
        labels=dict(genome="Number of Genes", index="Genome"),
        title="Genome Size Distribution"
    )
    fig.write_html("genome_frequency_bars.html", include_plotlyjs="cdn")


def alignment_qc_metrics(aln: pd.DataFrame):

    colorscale = [
        [0, 'rgb(255,255,255)'],
        [0.01, 'rgb(255,255,255)'],
        [0.1, 'rgb(12,51,131)'],
        [0.25, 'rgb(10,136,186)'],
        [0.5, 'rgb(242,211,56)'],
        [0.75, 'rgb(242,143,56)'],
        [1, 'rgb(217,30,30)']
    ]

    fig = px.density_heatmap(
        data_frame=aln,
        x="coverage",
        y="pident",
        labels=dict(
            coverage="Gene Length Detected (%)",
            pident="Amino Acid Identity (%)"
        ),
        marginal_x="histogram",
        marginal_y="histogram",
        title="Gene Detection QC Metrics",
        color_continuous_scale=colorscale
    )
    fig.write_html("alignment_qc_metrics.html", include_plotlyjs="cdn")
    

def alignment_heatmap(aln: pd.DataFrame):

    # Make a wide DataFrame
    wide_df = aln.pivot_table(
        index="genome",
        columns="sseqid",
        values="pident",
        aggfunc=max
    ).fillna(0)

    # Sort the rows and columns (default is euclidean distance)
    wide_df = sort(wide_df)
    wide_df = sort(wide_df.T).T

    fig = px.imshow(
        img=wide_df,
        aspect='auto',
        color_continuous_scale="Blues",
        labels=dict(sseqid="Gene", genome="Genome")
    )
    fig.write_html("alignment_heatmap.html", include_plotlyjs="cdn")


def gene_frequency_genome_size(aln: pd.DataFrame):

    # Make a wide DataFrame
    wide_df = (
        aln.pivot_table(
            index="sseqid",
            columns="genome",
            values="pident",
            aggfunc=max
        )
        .fillna(0)
        > 0
    )

    # Count up the number of genomes per gene
    vc = wide_df.sum(axis=1)

    # Count up the number of genes per genome at different thresholds
    plot_df = pd.DataFrame([
        (
            dict(
                genome=genome,
                count=count,
                threshold=threshold
            )
        )
        for threshold in range(1, vc.max())
        for genome, count in (
            wide_df
            .loc[vc >= threshold]
            .sum()
            .items()
        )
    ])

    fig = px.box(
        data_frame=plot_df,
        x="threshold",
        y="count",
        points="all",
        labels=dict(
            threshold="Minimum Number of Genomes per Gene",
            count="Number of Genes per Genome at Threshold"
        ),
        hover_data=["genome"]
    )
    fig.write_html("gene_frequency_genome_size.html", include_plotlyjs="cdn")


def sort(df: pd.DataFrame, method="average", metric="euclidean"):
    return df.reindex(
        index=df.index.values[
            hierarchy.leaves_list(
                hierarchy.linkage(
                    df,
                    method=method,
                    metric=metric
                )
            )
        ]
    )


if __name__ == "__main__":
    run()
