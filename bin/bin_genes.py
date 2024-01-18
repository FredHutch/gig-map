#!/usr/bin/env python3
import anndata as ad
import click
import logging
from typing import Dict, List, Union
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


class GeneData(ad.AnnData):

    @classmethod
    def from_alignment_csv(
        cls,
        aln_fp,
        gene_annot_fp,
        min_coverage=0.95,
        min_identity=0.95,
        min_genomes_per_gene=1
    ):
        """Set up an AnnData object with gene data."""

        # Read in the table
        logger.info(f"Reading in {aln_fp}")
        aln = pd.read_csv(aln_fp)

        # Make sure that we have the expected columns
        cls.validate_cnames(
            aln,
            ["genome", "sseqid", "coverage", "pident"],
            aln_fp
        )

        # Add a column for presence/absence
        aln = aln.assign(present=True)

        # Filter the alignments
        aln = cls.filter_aln(
            aln,
            min_coverage=min_coverage,
            min_identity=min_identity,
            min_genomes_per_gene=min_genomes_per_gene
        )

        # Make the wide tables
        wide_dfs: Dict[str, pd.DataFrame] = {
            values: aln.pivot_table(
                index="genome",
                columns="sseqid",
                values=values,
                aggfunc=max
            ).fillna(fillna)
            for values, fillna in [
                ("coverage", 0),
                ("pident", 0),
                ("present", False)
            ]
        }

        # Read in the gene annotations
        logger.info(f"Reading in {gene_annot_fp}")
        gene_annot = pd.read_csv(gene_annot_fp)
        cls.validate_cnames(gene_annot, ["gene_id", "combined_name"], gene_annot_fp)
        gene_annot.set_index("gene_id", inplace=True)

        ngenomes = aln["genome"].unique().shape[0]
        ngenes = aln["sseqid"].unique().shape[0]
        logger.info(f"Filtered to {ngenes:,} genes across {ngenomes:,} genomes")

        # Make sure that the indexes are the same
        var_ix = wide_dfs["present"].columns.values
        obs_ix = wide_dfs["present"].index.values

        return cls(
            X=wide_dfs["present"],
            layers={
                kw: (
                    wide_dfs
                    [kw]
                    .reindex(index=obs_ix, columns=var_ix)
                )
                for kw in ["coverage", "pident"]
            },
            var=gene_annot.reindex(index=var_ix)
        )

    @staticmethod
    def validate_cnames(df: pd.DataFrame, cnames: List[str], filename: str):
        """Make sure that a DataFrame contains a set of columns."""
        missing_cnames = [
            cname for cname in cnames
            if cname not in df.columns.values
        ]

        msg = f"Missing columns expected in {filename}: {', '.join(missing_cnames)}"
        assert len(missing_cnames) == 0, msg

    @classmethod
    def filter_aln(
        cls,
        aln: pd.DataFrame,
        min_coverage=0.95,
        min_identity=0.95,
        min_genomes_per_gene=1
    ) -> pd.DataFrame:
        """Filter the genome alignments."""

        logger.info(f"Total number of alignments: {aln.shape[0]:,}")
        aln = cls.filter_aln_query(aln, f"coverage >= {min_coverage}")
        aln = cls.filter_aln_query(aln, f"pident >= {min_identity}")

        # Filter by the number of genomes that each gene is found in
        return cls.filter_aln_ngenomes(aln, min_genomes_per_gene=min_genomes_per_gene)

    @staticmethod
    def filter_aln_query(aln: pd.DataFrame, query: str) -> pd.DataFrame:
        logger.info(f"Filtering to {query}")
        aln = aln.query(f"{query}")
        logger.info(f"Filtered number of alignments: {aln.shape[0]:,}")
        return aln

    @staticmethod
    def filter_aln_ngenomes(
        aln: pd.DataFrame,
        min_genomes_per_gene=1
    ) -> pd.DataFrame:
        logger.info(f"Keeping genes found in >= {min_genomes_per_gene} genomes")

        vc = aln.reindex(columns=['sseqid', 'genome'])['sseqid'].value_counts()
        tot = vc.shape[0]
        tot_counts = vc.sum()
        vc = vc.loc[(vc >= min_genomes_per_gene)]
        logger.info(f"{vc.shape[0]:,} / {tot:,}% genes meeting threshold")
        logger.info(f"{round(100 * vc.sum() / tot_counts, 1)} total gene content meeting threshold")
        return aln.loc[aln['sseqid'].isin(vc.index.values)]

    def bin_genes(
        self,
        method="average",
        metric="jaccard",
        max_dist_genes=0.05,
        min_bin_size=1
    ):
        """Combine genes into bins."""

        logger.info(f"Maximum {metric} distance for binning genes: {max_dist_genes}")
        logger.info(f"Minimum bin size: {min_bin_size:,}")

        bins = self.linkage_cluster(
            self.to_df().T,
            method=method,
            metric=metric,
            max_dist=max_dist_genes,
            min_size=min_bin_size,
            name="Gene Bin"
        )

        n_clust = bins.unique().shape[0]
        n_genes = bins.shape[0]
        logger.info(f"Formed {n_clust:,} bins with {n_genes:,} genes")

        # Assign to the var/bin slot
        self.var["bin"] = bins

    @staticmethod
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

    def summarize_bins(self):
        """Report summaries of bins that were formed."""

        # Total number of genes in each bin
        bin_counts = self.var["bin"].value_counts()

        # Make a bargraph
        fig = px.bar(
            data_frame=bin_counts.reset_index(),
            x="index",
            y="bin",
            labels=dict(bin="Number of Genes", index="Bin"),
            title="Gene Bin Size Distribution",
            log_y=True
        )
        fig.write_html("gene_bin_size_bars.html")

        # Write out a table of variable annotations,
        # which contain the binned genes
        (
            self.var
            .assign(
                bin_ix=(
                    self.var["bin"]
                    .apply(
                        lambda s: (
                            int(s.replace("Bin ", ""))
                            if isinstance(s, str) else
                            999
                        )
                    )
                ),
                n_genomes=self.to_df().sum()
            )
            .sort_values(
                by=["bin_ix", "n_genomes"],
                ascending=[True, False]
            )
            .drop(columns=["bin_ix"])
            .to_csv("gene_bins.csv")
        )

    def summarize_genomes(self):

        # Total number of genes in each bin
        bin_counts = self.var["bin"].value_counts()

        # Summarize the number of genes from each bin for each genome
        genome_content = (
            self
            .to_df()
            .T
            .groupby(self.var["bin"])
            .sum()
            .reset_index()
            .melt(id_vars="bin")
            .rename(
                columns=dict(
                    variable="genome",
                    value="n_genes_detected"
                )
            )
            .query("n_genes_detected > 0")
            .sort_values(
                by=["genome", "n_genes_detected"],
                ascending=[True, False]
            )
            .assign(
                prop_genes_detected=(
                    lambda d: d["n_genes_detected"] / d["bin"].apply(bin_counts.get)
                )
            )
        )
        genome_content.to_csv("genome_content.long.csv", index=None)

        for kw in ["n_genes_detected", "prop_genes_detected"]:
            self.obsm[kw] = (
                genome_content
                .pivot_table(
                    index="genome",
                    columns="bin",
                    values=kw
                )
                .reindex(self.obs_names)
                .fillna(0)
            )

        # Summarize the proportion of genome content which
        # is captured in any of the bins
        prop_in_bin = pd.DataFrame(dict(
            in_bin=self.obsm["n_genes_detected"].sum(axis=1),
            total=self.to_df().sum(axis=1)
        )).assign(
            prop=lambda d: d["in_bin"] / d["total"]
        ).sort_values(
            by=["prop", "in_bin"],
            ascending=False
        )
        logger.info("Proportion of genome content contained in bins:")
        for kw, val in prop_in_bin["prop"].describe().items():
            logger.info(f"{kw}: {val}")

        fig = make_subplots(rows=3, cols=1, shared_xaxes=True)
        for ix, (kw, xlabel) in enumerate(
            [
                ("prop", "Proportion Genes in Bins"),
                ("in_bin", "Number Genes in Bins"),
                ("total", "Total Genes Detected")
            ]
        ):
            fig.add_trace(
                go.Bar(
                    x=prop_in_bin.index.values,
                    y=prop_in_bin[kw].values,
                    hovertext=prop_in_bin[kw].values,
                ),
                row=ix+1,
                col=1
            )
            # Set the axis labels
            (
                fig
                ['layout']
                ['yaxis' if ix == 0 else f'yaxis{ix+1}']
                ['title']
            ) = xlabel

        prop_in_bin.to_csv("prop_genome_in_bins.csv")
        fig.write_html("prop_genome_in_bins.html")

    def group_genomes(
        self,
        max_dist_genomes=0.05,
        method="average",
        metric="braycurtis"
    ) -> pd.DataFrame:

        logger.info(f"Maximum {method} distance: {max_dist_genomes}")

        # Assign genomes into groups
        self.obs["group"] = self.linkage_cluster(
            self.obsm["n_genes_detected"],
            method=method,
            metric=metric,
            max_dist=max_dist_genomes,
            min_size=1,
            prefix="Group ",
            name="Genome Group"
        )
        self.obs.to_csv("genome_groups.csv")

    def make_heatmaps(self):

        # All genes, all genomes
        heatmap(
            (
                self
                .to_df()
                .reindex(
                    columns=self.var["bin"].dropna().index
                )
                .applymap(int)
            ),
            index_colors=self.obs["group"],
            index_colors_title="Genome Group",
            columns_colors=self.var["bin"],
            columns_colors_title="Gene Bin",
            filename="heatmap.unbinned_genes_all_genomes.html",
            metric="jaccard"
        )

        # Binned genes, all genomes
        heatmap(
            self.obsm["prop_genes_detected"],
            index_colors=self.obs["group"],
            index_colors_title="Genome Group",
            columns_bars=self.var["bin"].value_counts(),
            filename="heatmap.binned_genes_all_genomes.html",
            metric="braycurtis"
        )

        # Binned genes, binned genomes
        heatmap(
            (
                self.obsm["prop_genes_detected"]
                .groupby(self.obs["group"])
                .mean()
            ),
            index_bars=self.obs["group"].value_counts(),
            columns_bars=self.var["bin"].value_counts(),
            filename="heatmap.binned_genes_grouped_genomes.html",
            metric="braycurtis"
        )



class FigureBuilder:

    _panels: List[dict]
    _subplot_kwargs: dict

    def __init__(
        self,
        vertical_spacing=0,
        horizontal_spacing=0
    ):

        self._subplot_kwargs = dict(
            start_cell="bottom-left",
            rows=0,
            cols=0,
            shared_xaxes=True,
            shared_yaxes=True,
            column_widths=[],
            row_heights=[],
            vertical_spacing=vertical_spacing,
            horizontal_spacing=horizontal_spacing
        )

        self._panels = []

    def _add_dim(self, dim, trace, size, **kwargs):
        assert dim in ["row", "col"], dim
        size_kws = dict(
            col="column_widths",
            row="row_heights"
        )

        self._subplot_kwargs[f"{dim}s"] += 1
        self._subplot_kwargs[size_kws[dim]] = [1.0]
        for kw in ["row", "col"]:
            if self._subplot_kwargs[f"{kw}s"] == 0:
                self._subplot_kwargs[f"{kw}s"] = 1
                self._subplot_kwargs[size_kws[kw]] = [1.0]

        if self._subplot_kwargs[f"{dim}s"] > 1:
            self._subplot_kwargs[size_kws[dim]][0] -= size
            self._subplot_kwargs[size_kws[dim]].append(size)

        self._panels.append(dict(
            trace=trace,
            row=int(self._subplot_kwargs["rows"]) if dim == "row" else 1,
            col=int(self._subplot_kwargs["cols"]) if dim == "col" else 1,
            **kwargs
        ))

    def add_row(self, trace, height=0.05, **kwargs):
        self._add_dim("row", trace, height, **kwargs)

    def add_col(self, trace, width=0.05, **kwargs):
        self._add_dim("col", trace, width, **kwargs)

    def write_html(self, filename):
        logger.info(self._subplot_kwargs)
        fig = make_subplots(**self._subplot_kwargs)
        for panel in self._panels:
            fig.add_trace(
                panel["trace"],
                row=panel["row"],
                col=panel["col"]
            )
            if panel.get("log_x", False):
                fig.update_xaxes(
                    type="log",
                    row=panel["row"],
                    col=panel["col"]
                )
            if panel.get("log_y", False):
                fig.update_yaxes(
                    type="log",
                    row=panel["row"],
                    col=panel["col"]
                )

        fig.update_layout(
            plot_bgcolor='rgba(0, 0, 0, 0)',
            paper_bgcolor='rgba(0, 0, 0, 0)'
        )
        fig.write_html(filename)


def heatmap(
    X: pd.DataFrame,
    index_colors: Union[pd.Series, None] = None,
    index_colors_title: Union[str, None] = None,
    columns_colors: Union[pd.Series, None] = None,
    columns_colors_title: Union[str, None] = None,
    index_bars: Union[pd.Series, None] = None,
    columns_bars: Union[pd.Series, None] = None,
    filename="heatmap.html",
    method="average",
    metric="euclidean",
    colorscale="Blues"
):

    # Set up a FigureBuilder
    fb = FigureBuilder()

    # Get the sorted axes
    obs_ix = sort_index(X, method=method, metric=metric)
    var_ix = sort_index(X.T, method=method, metric=metric)

    # Set up the central heatmap
    fb.add_col(
        go.Heatmap(
            x=var_ix,
            y=obs_ix,
            z=X.reindex(index=obs_ix, columns=var_ix).values,
            colorscale=colorscale
        )
    )

    # Optionally add the column colors
    if columns_colors is not None:
        columns_colors = columns_colors.reindex(index=var_ix)
        fb.add_row(
            go.Heatmap(
                x=var_ix,
                y=[
                    (
                        columns_colors.name
                        if columns_colors_title is None else
                        columns_colors_title
                    )
                ],
                z=alternate_colors(columns_colors).to_frame().T.values,
                text=columns_colors.to_frame().T.values,
                colorscale="bluered",
                showscale=False
            )
        )

    # Optionally add the row colors
    if index_colors is not None:
        index_colors = index_colors.reindex(index=obs_ix)
        fb.add_col(
            go.Heatmap(
                y=obs_ix,
                x=[
                    (
                        index_colors.name
                        if index_colors_title is None else
                        index_colors_title
                    )
                ],
                z=alternate_colors(index_colors).to_frame().values,
                text=index_colors.to_frame().values,
                colorscale="bluered",
                showscale=False
            )
        )

    # Optionally add marginal bargraphs
    if columns_bars is not None:
        columns_bars = columns_bars.reindex(index=var_ix)

        fb.add_row(
            go.Bar(
                x=var_ix,
                y=columns_bars.values,
                hovertext=columns_bars.values,
                showlegend=False
            ),
            log_y=True
        )

    if index_bars is not None:
        index_bars = index_bars.reindex(index=obs_ix)

        fb.add_col(
            go.Bar(
                y=obs_ix,
                x=index_bars.values,
                hovertext=index_bars.values,
                orientation="h",
                showlegend=False
            )
        )

    fb.write_html(filename)


def alternate_colors(vals: pd.Series, n=2) -> pd.Series:
    output = []
    for ix, val in enumerate(vals.values):
        if ix == 0:
            output.append(0)
        elif val != vals.values[ix-1]:
            output.append((output[-1] + 1) % n)
        else:
            output.append(output[-1])

    return pd.Series(output, index=vals.index.values)


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


@click.command
@click.option('--genome_aln', type=str)
@click.option('--gene_annot', type=str)
@click.option('--min_coverage', type=float)
@click.option('--min_identity', type=float)
@click.option('--min_genomes_per_gene', type=int)
@click.option('--max_dist_genes', type=float)
@click.option('--min_bin_size', type=int)
@click.option('--max_dist_genomes', type=float)
def main(
    genome_aln,
    gene_annot,
    min_coverage,
    min_identity,
    min_genomes_per_gene,
    max_dist_genes,
    min_bin_size,
    max_dist_genomes
):

    # Read in the genome alignment data and set up an object
    gene_data = GeneData.from_alignment_csv(
        genome_aln,
        gene_annot,
        min_coverage=min_coverage,
        min_identity=min_identity,
        min_genomes_per_gene=min_genomes_per_gene
    )

    # bin the genes
    gene_data.bin_genes(
        max_dist_genes=max_dist_genes,
        min_bin_size=min_bin_size
    )

    # Summarize the bins
    gene_data.summarize_bins()

    # Summarize the genomes
    gene_data.summarize_genomes()

    # Group the genomes
    gene_data.group_genomes(
        max_dist_genomes=max_dist_genomes
    )

    # Make heatmaps
    gene_data.make_heatmaps()


if __name__ == "__main__":
    main()
