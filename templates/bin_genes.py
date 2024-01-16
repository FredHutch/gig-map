#!/usr/bin/env python3
from collections import defaultdict
import logging
from typing import Dict, List
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.cluster import hierarchy

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [bin_genes.py] %(message)s'
)
logger = logging.getLogger('bin_genes.py')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Write to file
fileHandler = logging.FileHandler("bin_genes.log")
fileHandler.setFormatter(logFormatter)
fileHandler.setLevel(logging.INFO)
logger.addHandler(fileHandler)


def validate_cnames(df: pd.DataFrame, cnames: List[str], filename: str):
    """Make sure that a DataFrame contains a set of columns."""
    missing_cnames = [
        cname for cname in cnames
        if cname not in df.columns.values
    ]

    msg = f"Missing columns expected in {filename}: {', '.join(missing_cnames)}"
    assert len(missing_cnames) == 0, msg


def main():

    # Read in the genome alignments
    logger.info("Reading in ${genome_aln}")
    aln = pd.read_csv("${genome_aln}")
    validate_cnames(
        aln,
        ["genome", "sseqid", "coverage", "pident"],
        "${genome_aln}"
    )

    # Read in the gene annotations
    logger.info("Reading in ${gene_annot}")
    gene_annot = pd.read_csv("${gene_annot}")
    validate_cnames(gene_annot, ["gene_id", "combined_name"], "${gene_annot}")

    # Filter the alignments
    aln = filter_aln(aln)

    ngenomes = aln["genome"].unique().shape[0]
    ngenes = aln["sseqid"].unique().shape[0]
    logger.info(f"Filtered to {ngenes:,} genes across {ngenomes:,} genomes")

    # bin the genes
    bins = bin_genes(aln)

    # Summarize the bins
    summarize_bins(bins, gene_annot)

    # Summarize the genomes
    genome_content = summarize_genomes(bins, aln)

    # Group the genomes
    genome_groups = group_genomes(genome_content)

    # Make heatmaps
    make_heatmaps(aln, bins, genome_content, genome_groups)


def filter_aln(aln: pd.DataFrame) -> pd.DataFrame:
    """Filter the genome alignments."""

    logger.info(f"Total number of alignments: {aln.shape[0]:,}")
    aln = filter_aln_query(aln, "coverage >= ${params.min_coverage}")
    aln = filter_aln_query(aln, "pident >= ${params.min_identity}")

    # Filter by the number of genomes that each gene is found in
    return filter_aln_ngenomes(aln)


def filter_aln_query(aln: pd.DataFrame, query: str) -> pd.DataFrame:
    logger.info(f"Filtering to {query}")
    aln = aln.query(f"{query}")
    logger.info(f"Filtered number of alignments: {aln.shape[0]:,}")
    return aln


def filter_aln_ngenomes(aln: pd.DataFrame) -> pd.DataFrame:
    min_genomes_per_gene = int("${params.min_genomes_per_gene}")
    logger.info(f"Keeping genes found in >= {min_genomes_per_gene} genomes")

    vc = aln.reindex(columns=['sseqid', 'genome'])['sseqid'].value_counts()
    tot = vc.shape[0]
    tot_counts = vc.sum()
    vc = vc.loc[(vc >= min_genomes_per_gene)]
    logger.info(f"{vc.shape[0]:,} / {tot:,}% genes meeting threshold")
    logger.info(f"{round(100 * vc.sum() / tot_counts, 1)} total gene content meeting threshold")
    return aln.loc[aln['sseqid'].isin(vc.index.values)]


def bin_genes(
    aln: pd.DataFrame,
    method="average",
    metric="jaccard"
) -> pd.Series:
    """Combine genes into bins."""

    # Make a wide DataFrame
    wide_df = (
        aln
        .assign(present=True)
        .reindex(
            columns=["genome", "sseqid", "present"]
        )
        .drop_duplicates()
    ).pivot(
        columns="genome",
        index="sseqid",
        values="present"
    ).fillna(
        False
    )

    max_dist_genes = float("${params.max_dist_genes}")
    logger.info(f"Maximum {metric} distance for binning genes: {max_dist_genes}")

    min_bin_size = int("${params.min_bin_size}")
    logger.info(f"Minimum bin size: {min_bin_size:,}")

    bins = linkage_cluster(
        wide_df,
        method=method,
        metric=metric,
        max_dist=max_dist_genes,
        min_size=min_bin_size,
        name="Gene Bin"
    )

    n_clust = bins.unique().shape[0]
    n_genes = bins.shape[0]
    logger.info(f"Formed {n_clust:,} bins with {n_genes:,} genes")

    return bins


def linkage_cluster(
    wide_df: pd.DataFrame,
    method="average",
    metric="euclidean",
    max_dist=0.05,
    min_size=1,
    prefix="Bin ",
    name="bin"
) -> pd.Series:

    intro = "Performing linkage cluster analysis: "
    kwarg_str = ", ".join([
        f"{kw}={val}"
        for kw, val in dict(
            method=method,
            metric=metric,
            max_dist=max_dist,
        ).items()
    ])
    logger.info(f"{intro}{wide_df.shape[0]:,} items ({kwarg_str})")
    L = hierarchy.linkage(
        wide_df,
        method=method,
        metric=metric
    )
    bins = hierarchy.fcluster(
        L,
        max_dist,
        criterion="distance"
    )

    # Count the number of genes in each bin
    vc = pd.Series(bins).value_counts().sort_values(ascending=False)
    logger.info(f"Formed {vc.shape[0]:,} bins in total (before filtering)")

    # Rename the bins, and mask anything below the threshold
    bin_map = {
        old_ix: f"{prefix}{new_ix + 1}"
        for new_ix, (old_ix, count) in enumerate(vc.items())
        if count >= min_size
    }

    # Map each gene name to the bin assignment
    bins = (
        pd.DataFrame([
            {
                "index": index,
                name: bin_map.get(bin_ix)
            }
            for index, bin_ix in zip(
                wide_df.index.values,
                bins
            )
        ])
        .dropna()
        .set_index("index")
        [name]
    )

    return bins


def summarize_bins(
    bins: pd.Series,
    gene_annot: pd.DataFrame
):
    """Report summaries of bins that were formed."""

    # Total number of genes in each bin
    bin_counts = bins.value_counts()

    # Make a bargraph
    fig = px.bar(
        data_frame=bin_counts.reset_index(),
        x="index",
        y="Gene Bin",
        labels=dict(bin="Number of Genes", index="Bin"),
        title="Gene Bin Size Distribution",
        log_y=True
    )
    fig.write_html("gene_bin_size_bars.html")

    # Make the bins into a table
    bins = pd.DataFrame(dict(bin=bins))

    # Add the annotations (if any)
    if gene_annot.shape[0] > 0:
        logger.info("Adding gene annotations")
        bins = bins.merge(
            gene_annot,
            left_index=True,
            right_on="gene_id",
            how="left"
        )
        n = bins['combined_name'].dropna().shape[0]
        logger.info(f"Genes in bins with annotations: {n:,}")

    # Write out a table of the binned genes
    bins.to_csv("bins.csv.gz")


def summarize_genomes(bins: pd.DataFrame, aln: pd.DataFrame):

    # Total number of genes in each bin
    bin_counts = bins.value_counts()

    # Add the bin annotation to the alignments table
    aln = aln.assign(
        bin=(
            aln
            ['sseqid']
            .apply(
                bins.get
            )
        )
    )

    # Summarize the number of genes from each bin for each genome
    genome_content = pd.DataFrame([
        dict(
            genome=genome,
            bin=bin,
            n_genes=df.shape[0]
        )
        for (genome, bin), df in (
            aln
            .reindex(
                columns=['genome', 'bin', 'sseqid']
            )
            .drop_duplicates()
            .groupby(['genome', 'bin'])
        )
    ]).sort_values(
        by=["genome", "n_genes"],
        ascending=[True, False]
    )

    # Calculate the proportion of each bin which is found
    genome_content = genome_content.assign(
        bin_prop=genome_content.apply(
            lambda r: r['n_genes'] / bin_counts[r["bin"]],
            axis=1
        )
    )
    genome_content.to_csv("genome_content.long.csv", index=None)

    return genome_content


def group_genomes(
    genome_content: pd.DataFrame,
    method="average",
    metric="braycurtis"
) -> pd.DataFrame:

    # Make a wide DataFrame
    wide = {
        values: (
            genome_content
            .pivot_table(
                index="genome",
                columns="bin",
                values=values
            )
            .fillna(0)
        )
        for values in ["bin_prop", "n_genes"]
    }

    # Order the genomes
    genome_order = sort_index(
        wide["n_genes"],
        method=method,
        metric=metric
    )
    for kw in wide:
        wide[kw] = wide[kw].reindex(index=genome_order)

    max_dist_genomes = float("${params.max_dist_genomes}")
    logger.info(
        f"Maximum {metric} distance for grouping genomes: {max_dist_genomes}"
    )

    # Assign genomes into groups
    genome_groups = linkage_cluster(
        wide["n_genes"],
        method=method,
        metric=metric,
        max_dist=max_dist_genomes,
        min_size=1,
        prefix="Group ",
        name="Genome Group"
    )

    return genome_groups


def sort_index(df: pd.DataFrame, method="average", metric="euclidean"):

    intro = "Sorting index: "
    kwarg_str = ", ".join([
        f"{kw}={val}"
        for kw, val in dict(
            method=method,
            metric=metric,
        ).items()
    ])
    logger.info(f"{intro}{df.shape[0]:,} items ({kwarg_str})")

    return df.index.values[
        hierarchy.leaves_list(
            hierarchy.linkage(
                df,
                method=method,
                metric=metric
            )
        )
    ]


def make_heatmaps(
    aln: pd.DataFrame,
    bins: pd.DataFrame,
    genome_content: pd.DataFrame,
    genome_groups: pd.DataFrame
):
    # All genes, all genomes
    heatmap(
        aln.assign(present=1),
        index="genome",
        index_colors=genome_groups,
        columns="sseqid",
        columns_colors=bins,
        values="present",
        filename="heatmap.unbinned_genes_all_genomes.html",
        metric="jaccard"
    )

    # Binned genes, all genomes
    heatmap(
        genome_content,
        index="genome",
        index_colors=genome_groups,
        columns="bin",
        columns_colors=None,
        values="bin_prop",
        filename="heatmap.binned_genes_all_genomes.html",
        metric="braycurtis"
    )

    # Binned genes, binned genomes
    heatmap(
        (
            genome_content
            .assign(
                genome_group=genome_content["genome"].apply(
                    genome_groups.get
                )
            )
            .groupby(
                ["genome_group", "bin"]
            )
            ["bin_prop"]
            .mean()
            .reset_index()
        ),
        index="genome_group",
        index_colors=None,
        columns="bin",
        columns_colors=None,
        values="bin_prop",
        filename="heatmap.binned_genes_binned_genomes.html",
        metric="braycurtis"
    )


def heatmap(
    long_df: pd.DataFrame,
    index: str,
    index_colors: pd.Series,
    columns: str,
    columns_colors: pd.Series,
    values: str,
    filename: str,
    fillna=0,
    aggfunc=max,
    method="average",
    metric="euclidean",
    colorscale="Blues"
):
    logger.info(f"Preparing heatmap for {filename}")
    for kw, cname in [
        ('index', index),
        ('columns', columns),
        ('values', values)
    ]:
        opts = ', '.join(long_df.columns.values)
        msg = f"Could not find {kw} column {cname} in {opts}"
        assert cname in long_df.columns.values, msg

    # Make a wide table
    wide_df = long_df.pivot_table(
        index=index,
        columns=columns,
        values=values,
        aggfunc=aggfunc
    ).fillna(fillna)

    nrows, ncols = wide_df.shape
    logger.info(
        f"Prepared wide table with {nrows:,} rows and {ncols:,} columns"
    )

    if index_colors is not None:
        # Filter to rows/cols with values for annotations
        index_overlap = list(
            set(wide_df.index.values) &
            set(index_colors.index.values)
        )
        msg = f"No shared values found between {index} and the assigned colors"
        assert len(index_overlap) > 0, msg

        wide_df = wide_df.reindex(
            index=index_overlap
        )
        index_colors = pd.DataFrame({index_colors.name: index_colors})

    if columns_colors is not None:
        columns_overlap = list(
            set(wide_df.columns.values) &
            set(columns_colors.index.values)
        )
        msg = f"No shared values found between {columns} and the assigned colors"
        assert len(columns_overlap) > 0, msg

        wide_df = wide_df.reindex(
            columns=columns_overlap
        )

        columns_colors = pd.DataFrame({columns_colors.name: columns_colors})

    # Set up subplot kwargs
    # Default setup is only really displaying a single plot
    subplot_kwargs = dict(
        rows=2,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        column_widths=[1, 0],
        row_heights=[1, 0],
        vertical_spacing=0,
        horizontal_spacing=0
    )

    # Sort both axes
    wide_df = wide_df.reindex(
        index=sort_index(wide_df, method=method, metric=metric),
        columns=sort_index(wide_df.T, method=method, metric=metric)
    )
    if index_colors is not None:
        index_colors = index_colors.reindex(index=wide_df.index.values)
        subplot_kwargs["column_widths"] = [0.9, 0.1]
        subplot_kwargs["vertical_spacing"] = 0.02

    if columns_colors is not None:
        columns_colors = columns_colors.reindex(index=wide_df.columns.values)
        subplot_kwargs["row_heights"] = [0.1, 0.9]
        subplot_kwargs["horizontal_spacing"] = 0.02

    # Set up the subplots
    fig = make_subplots(**subplot_kwargs)

    # Central heatmap
    fig.add_trace(
        go.Heatmap(
            x=wide_df.columns.values,
            y=wide_df.index.values,
            z=wide_df.values,
            colorscale=colorscale
        ),
        row=2,
        col=1
    )

    if columns_colors is not None:

        # Column colors
        fig.add_trace(
            go.Heatmap(
                x=columns_colors.index.values,
                y=columns_colors.columns.values,
                z=(
                    alternate_colors(columns_colors)
                    .T
                    .values
                ),
                text=columns_colors.T.values,
                colorscale="bluered",
                showscale=False
            ),
            row=1,
            col=1
        )

    if index_colors is not None:

        # Index colors
        fig.add_trace(
            go.Heatmap(
                x=index_colors.columns.values,
                y=index_colors.index.values,
                z=(
                    alternate_colors(index_colors)
                    .values
                ),
                text=index_colors.values,
                colorscale="bluered",
                showscale=False
            ),
            row=2,
            col=2
        )

    fig.write_html(filename)


def alternate_colors(vals: pd.DataFrame, n=2) -> pd.DataFrame:
    output: Dict[str, List[int]] = {}
    for cname, cvals in vals.items():
        output[cname] = []
        for ix, val in enumerate(cvals.values):
            if ix == 0:
                output[cname].append(0)
            elif val != cvals.values[ix-1]:
                output[cname].append((output[cname][-1] + 1) % n)
            else:
                output[cname].append(output[cname][-1])

    return pd.DataFrame(output, index=vals.index.values)


if __name__ == "__main__":
    main()
