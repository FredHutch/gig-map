#!/usr/bin/env python3
import gzip
import json
from pathlib import Path
from typing import Tuple
import click
import numpy as np
import pandas as pd
import sys
import logging
from scipy.cluster import hierarchy
from plotly import graph_objects as go
from multiprocessing import Pool
logger = logging.getLogger("plot_bins.py")
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)


@click.command
@click.option('--genome_aln', type=click.Path(exists=True))
@click.option('--gene_bins', type=click.Path(exists=True))
@click.option('--min_coverage', type=float)
@click.option('--min_identity', type=float)
@click.option('--threads', type=int, default=1)
def main(genome_aln, gene_bins, min_coverage, min_identity, threads):
    # Read in the table of alignments
    logger.info("Reading in the alignments table: %s", genome_aln)
    aln = pd.read_csv(genome_aln, index_col=0)
    logger.info("Done reading in the alignments table: %s", aln.shape)

    # Filter on minimum coverage and identity
    logger.info(f"Filtering on minimum coverage: {min_coverage}")
    aln = aln.query(f"coverage >= {min_coverage}")
    logger.info(f"Alignments passing filter: {aln.shape[0]:,}")
    logger.info(f"Filtering on minimum identity: {min_identity}")
    aln = aln.query(f"pident >= {min_identity}")
    logger.info(f"Alignments passing filter: {aln.shape[0]:,}")

    # Read in the assignment of genes to bins
    bins = pd.read_csv(gene_bins)
    logger.info("Read in the bins table: %s", bins.shape)

    # Drop any genes that have NaN in the bin column
    bins = bins.dropna(subset=["bin"])

    logger.info(f"Number of genes in bins: {bins.shape[0]:,}")
    logger.info(f'Number of bins: {bins["bin"].nunique():,}')

    # Merge the two tables
    merged = aln.merge(
        bins,
        how="inner",
        left_on="sseqid",
        right_on="gene_id"
    ).rename(
        columns=dict(
            qseqid="contig",
            qstart="contig_start",
            qend="contig_end",
            sstart="gene_start",
            send="gene_end",
            combined_name="gene_name"
        )
    ).reindex(columns=[
        "genome",
        "contig",
        "contig_start",
        "contig_end",
        "gene_id",
        "gene_start",
        "gene_end",
        "pident",
        "coverage",
        "bin",
        "gene_name"
    ])

    # Set up a multiprocessing pool
    pool = Pool(threads)

    # For each bin, make a display
    pool.map(
        layout_bin,
        [
            (bin, aln_df)
            for bin, aln_df in merged.groupby("bin")
        ]
    )


class Layout:
    """Object used to coordinate the layout of genes in a bin across genomes"""
    # Keep track of which genes have been added for which genomes
    # items: tuple[genome, gene]
    added_genes: set

    # Build a DataFrame with the global coordnates for each gene
    # in the layout
    coords: pd.DataFrame

    # Working copy of the alignments that need to be added
    aln_df: pd.DataFrame

    def __init__(self, aln_df: pd.DataFrame):
        self.added_genes = set()
        self.coords = pd.DataFrame()
        self.aln_df = aln_df

    def get_next_contig_to_add(self) -> Tuple[str, str]:
        """Return the next genome/contig to add to the layout."""

        # Compute the number of genes that have not been added
        # which are in each contig
        contig_priority = self.find_contig_priority()

        # If there are no genes that have not been added
        # then we are done
        if (
            contig_priority is None or
            contig_priority.shape[0] == 0 or
            contig_priority["n_pending"].max() == 0
        ):
            return (None, None)

        top_genome = contig_priority.iloc[0]['genome']
        top_contig = contig_priority.iloc[0]['contig']

        return (top_genome, top_contig)

    def add_all_contigs(self):
        """
        If there is a contig in the alignments which has
        genes that have not yet been added, find the contig
        which has the most genes that have not been added
        and add that to the global layout.
        """
        top_genome, top_contig = self.get_next_contig_to_add()

        while top_genome is not None:
            # Add the top contig to the layout
            self.add_contig(top_genome, top_contig)

            # Drop that contig from the running alignments
            self.aln_df = (
                self.aln_df
                .query(f"genome != '{top_genome}' or contig != '{top_contig}'")
            )
            # Get the next set
            top_genome, top_contig = self.get_next_contig_to_add()

    def add_contig(self, top_genome: str, top_contig: str):
        """Add a single contig to the layout."""

        contig_df = (
            self.aln_df
            .query(f"genome == '{top_genome}'")
            .query(f"contig == '{top_contig}'")
        )
        logger.info(f"Adding contig '{top_contig}' from genome '{top_genome}'")

        # Generate the potential coordinates for each gene
        # in the contig using the forward and reverse strand,
        # and the global layout coordinates
        fwd_coords, fwd_error = self.generate_coords(contig_df, "fwd")
        rev_coords, rev_error = self.generate_coords(contig_df, "rev")

        # Use the coordinates with the lower error
        if fwd_error <= rev_error:
            self.coords = pd.concat([self.coords, fwd_coords])
        else:
            self.coords = pd.concat([self.coords, rev_coords])

        # Record which genes have been added
        self.added_genes.update(set(zip(contig_df["genome"], contig_df["gene_id"])))

    def generate_coords(self, contig_df: pd.DataFrame, orientation: str):
        """
        Compute the coordinates of the contig in the global layout.
        """
        # If the orientation is rev, invert the coordinates
        if orientation == 'rev':
            contig_df = contig_df.assign(
                contig_start=contig_df["contig_start"] * -1,
                contig_end=contig_df["contig_end"] * -1,
            )

        # Subtract the mean offset to the pre-existing median coordinates
        # for each gene.
        # Also compute the difference between that adjusted coordinate
        # and the median coordinate for each gene
        global_contig_df = self.add_global_pos(contig_df)

        # Compute the mean absolute error
        if global_contig_df["error_start"].isnull().all():
            mean_error = 0
        else:
            mean_error = np.mean([
                global_contig_df["error_start"].dropna().abs().mean(),
                global_contig_df["error_end"].dropna().abs().mean(),
            ])

        # Clean up the global coordinates before returning them,
        # since that table will be added to self.coords directly
        global_contig_df = (
            global_contig_df
            .drop(columns=[
                "error_start",
                "error_end",
            ])
            .assign(
                contig_start=global_contig_df["contig_start"].abs(),
                contig_end=global_contig_df["contig_end"].abs(),
            )
        )
        if "offset_start" in global_contig_df.columns:
            global_contig_df = global_contig_df.drop(
                columns=[
                    "offset_start",
                    "offset_end",
                    "median_global_start",
                    "median_global_end",
                ]
            )

        return global_contig_df, mean_error

    def add_global_pos(self, contig_df: pd.DataFrame):
        """
        Compute the median offset for every gene to the existing coordinates.
        Also compute the difference between the adjusted position and that median.
        """

        # If no global coordinates are present
        if self.coords.shape[0] == 0:
            # Then just start at the beginning
            offset_val = min(
                contig_df["contig_start"].min(),
                contig_df["contig_end"].min(),
            )
            # Move the coordinates to start at the beginning
            # And there is no error to consider
            global_df = contig_df.assign(
                global_start=(contig_df["contig_start"] - offset_val).apply(int),
                global_end=(contig_df["contig_end"] - offset_val).apply(int),
                error_start=0,
                error_end=0
            )
        else:
            # Get the median positions for each gene
            median_positions = self.median_positions()

            # Compute the offset for every gene in the contig
            offset_df = (
                contig_df
                .merge(median_positions, on="gene_id", how="outer")
                .assign(
                    offset_start=lambda d: d["contig_start"] - d["median_global_start"],
                    offset_end=lambda d: d["contig_end"] - d["median_global_end"],
                )
            )

            # If there are no shared positions to use for the overlap
            if offset_df["offset_start"].isnull().all():
                # Then add this new contig to the end
                offset_val = min(
                    contig_df["contig_start"].min(),
                    contig_df["contig_end"].min(),
                ) + self.global_end() + 100

            # Use the median
            offset_val = np.mean([
                offset_df["offset_start"].dropna().median(),
                offset_df["offset_end"].dropna().median(),
            ])

            # Adjust the coordinates and add the difference between
            # the median and this new adjusted position
            global_df = offset_df.assign(
                global_start=offset_df["contig_start"] - offset_val,
                global_end=offset_df["contig_end"] - offset_val,
                error_start=lambda d: d["median_global_start"] - d["global_start"],
                error_end=lambda d: d["median_global_end"] - d["global_end"],
            )

        return global_df

    # Compute median position for each gene
    # in the global layout coordinates
    def median_positions(self) -> pd.DataFrame:
        return (
            self.coords
            .reindex(columns=["gene_id", "global_start", "global_end"])
            .groupby("gene_id")
            .median()
            .rename(columns=dict(
                global_start="median_global_start",
                global_end="median_global_end",
            ))
        )

    def global_end(self) -> int:
        """Find the last position used in global space."""
        return max(
            self.coords["global_start"].max(),
            self.coords["global_end"].max(),
        )

    def find_contig_priority(self) -> pd.DataFrame:
        """
        Compute the number of genes that have not been added
        which are in each contig
        """
        if self.aln_df.shape[0] == 0:
            return

        return pd.DataFrame([
            dict(
                genome=genome,
                contig=contig,
                n_pending=len(set([(genome, gene) for gene in df["gene_id"]]) - self.added_genes),
                n_total=df["gene_id"].nunique()
            )
            for (genome, contig), df in self.aln_df.groupby(["genome", "contig"])
        ]).sort_values(
            by=["n_pending", "n_total"],
            ascending=False
        )

    def save_csv(self, bin: str):
        """Save the coordinates to a CSV"""

        output_fp = Path(".") / "layout" / f"{bin}.csv.gz"
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving bin coordinates to {output_fp}")
        self.coords.to_csv(output_fp, index=None)

    def plot(self, bin: str):
        logger.info(f"Generating display for {bin}")
        genome_order = self.sort_genomes()

        # Set up a figure
        # fig = make_subplots(specs=[[{"secondary_x": True}]])
        fig = go.Figure()

        for genome_ix, genome_name in enumerate(genome_order):
            for contig_name, contig_df in (
                self.coords
                .query(f"genome == '{genome_name}'")
                .groupby('contig')
            ):
                self.plot_contig(
                    genome_ix,
                    genome_name,
                    contig_name,
                    contig_df,
                    fig
                )

        # add a second set of labels across the top with the gene names
        gene_labels = (
            self.coords
            .groupby(["gene_id", "gene_name"])
            .apply(lambda d: np.mean([d["global_start"].mean(), d["global_end"].mean()]), include_groups=False)
            .reset_index()
            .rename(columns={0: "global_middle"})
            .assign(
                text=lambda d: d.apply(
                    lambda r: f'{r["gene_id"]} {r["gene_name"]}',
                    axis=1
                )
            )
        )

        fig.update_layout(
            template="simple_white",
            xaxis=dict(
                automargin=True,
                tickvals=gene_labels["global_middle"].tolist(),
                ticktext=gene_labels["text"].tolist(),
                tickmode="array",
            ),
            yaxis=dict(
                tickvals=list(range(genome_order.shape[0])),
                ticktext=list(genome_order),
                tickmode="array",
            )
        )

        # Write as HTML
        output_fp = Path(".") / "layout" / f"{bin}.html"
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving bin layout image to {output_fp}")
        fig.write_html(output_fp)
        assert output_fp.exists()

        # Write as JSON
        output_fp = Path(".") / "layout" / f"{bin}.json.gz"
        logger.info(f"Saving bin layout image JSON to {output_fp}")
        with gzip.open(output_fp, 'wt') as handle:
            handle.write(json.dumps(fig.to_plotly_json()))
        assert output_fp.exists()

    def plot_contig(
        self,
        genome_ix: int,
        genome_name: str,
        contig_name: str,
        contig_df: pd.DataFrame,
        fig: go.Figure
    ):

        # Draw each gene
        for _, r in contig_df.iterrows():
            # If we are in the forward direction, draw above the line, otherwise below
            fwd = (r["global_start"] < r["global_end"])
            height = 0.2 * (1 if fwd else -1)

            hovertext = "<br>".join([
                f"{kw.title().replace('_', ' ')}: {r[kw]}"
                for kw in ["gene_id", "gene_name", "genome", "contig", "contig_start", "contig_end"]
            ])

            fig.add_trace(
                go.Scatter(
                    mode="lines",
                    x=[r["global_start"], r["global_end"], r["global_end"], r["global_start"]],
                    y=[genome_ix, genome_ix, genome_ix + height, genome_ix + height],
                    fillcolor="blue" if fwd else "orange",
                    showlegend=False,
                    fill="toself",
                    marker_color="blue" if fwd else "orange",
                    marker_line_width=0
                )
            )
            # Add a dummy point for the hovertext
            fig.add_trace(
                go.Scatter(
                    x=[np.mean([r["global_start"], r["global_end"]])],
                    y=[genome_ix + (height / 2)],
                    text=[hovertext],
                    hovertemplate="%{text}<extra></extra>",
                    mode="none",
                    showlegend=False
                )
            )

        # Draw a line for the contig
        contig_min = min(contig_df["global_start"].min(), contig_df["global_end"].min())
        contig_max = max(contig_df["global_start"].max(), contig_df["global_end"].max())
        fig.add_trace(
            go.Scatter(
                x=[contig_min, contig_max],
                y=[genome_ix, genome_ix],
                name=contig_name,
                showlegend=False,
                mode="lines",
                text=f"{genome_name} ({contig_name})",
                hovertemplate="%{text}<extra></extra>",
                marker=dict(
                    color="grey"
                )
            )
        )

    def sort_genomes(self):
        """Find the order of genomes with linkage clustering."""

        if self.coords["genome"].nunique() < 3:
            return self.coords["genome"].unique()

        # Make a wide table with the gene start coordinates
        wide_df = self.coords.pivot_table(
            index="genome",
            columns="gene_id",
            values="global_start",
            aggfunc="min"
        ).fillna(0)

        # Return the sorted list
        return wide_df.index.values[
            hierarchy.leaves_list(
                hierarchy.linkage(
                    wide_df.values,
                    method="ward",
                    metric="euclidean"
                )
            )
        ]


def layout_bin(inputs: Tuple[str, pd.DataFrame]):
    """Generate the layout for a bin across every genome that it is found in"""

    # Split up the inputs
    bin, aln_df = inputs

    # If there is only one gene or more than 500, skip it
    if (
        aln_df["gene_id"].nunique() == 1 or
        aln_df["gene_id"].nunique() > 500
    ):
        return

    logger.info(f"Computing layout for {bin}")

    # Make the layout object
    layout = Layout(aln_df)

    # Add the alignment information
    layout.add_all_contigs()
    logger.info(''.join([
        f"Found positions for {aln_df['genome'].nunique():,} genomes, "
        f"{aln_df['contig'].nunique():,} contigs, and "
        f"{aln_df['gene_id'].nunique():,} genes."
    ]))

    # Save the table
    layout.save_csv(bin)

    # Plot the coordinates
    layout.plot(bin)


if __name__ == "__main__":
    main()
