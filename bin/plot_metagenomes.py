#!/usr/bin/env python3
from typing import List
import anndata as ad
import click
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.cluster import hierarchy
from statsmodels.stats.multitest import multipletests

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [plot_metagenomes.py] %(message)s'
)
logger = logging.getLogger('plot_metagenomes.py')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Write to file
fileHandler = logging.FileHandler("plot_metagenomes.log")
fileHandler.setFormatter(logFormatter)
fileHandler.setLevel(logging.INFO)
logger.addHandler(fileHandler)


def fix_name(n: str, options: List[str]):
    if n in options:
        return n
    for m in options:
        if m.replace(" ", ".") == n:
            return m
    return n


class Metagenome:

    adata: ad.AnnData

    def __init__(
        self,
        stats,
        h5ad
    ):

        stats = self._read_csv(stats).set_index("feature")
        self.adata = ad.read_h5ad(h5ad)

        # Restore the variable names that may have been mangled by R
        stats = stats.rename(
            index=lambda n: fix_name(n, self.adata.var_names)
        )
        self.adata.var = self.adata.var.merge(
            stats,
            left_index=True,
            right_index=True,
            how="outer"
        )
        self.adata.var = self.adata.var.assign(**{
            kw: self.adata.var[kw].fillna(val)
            for kw, val in [("estimate", 0), ("std_error", 0), ("p_value", 0.99)]
        })
        self.adata.var["qvalue"] = multipletests(self.adata.var["p_value"].fillna(1), 0.1, "fdr_bh")[1]
        self.adata.var["neg_log10_pvalue"] = -np.log10(self.adata.var["p_value"])
        self.adata.var["neg_log10_qvalue"] = -np.log10(self.adata.var["qvalue"])
        for line in str(self.adata).split("\n"):
            logger.info(line)

    def _read_csv(self, fp, **kwargs):
        logger.info(f"Reading in {fp}")
        df: pd.DataFrame = pd.read_csv(fp, **kwargs)
        logger.info(f"Read in {df.shape[0]:,} rows and {df.shape[1]:,} columns")
        self.log_df(df.head().T.head().T)
        return df

    @staticmethod
    def log_df(df: pd.DataFrame, **kwargs):
        for line in df.to_csv(**kwargs).split("\n"):
            logger.info(line)

    @staticmethod
    def log_scale(df: pd.DataFrame):

        lowest = df.apply(lambda c: c[c > 0].min()).min()
        return df.clip(lower=lowest).apply(np.log10)

    def sort_index(self, df: pd.DataFrame, metric="cosine", method="average"):
        try:
            return df.index.values[
                hierarchy.leaves_list(
                    hierarchy.linkage(
                        df.values,
                        metric=metric,
                        method=method
                    )
                )
            ]
        except ValueError as e:
            if metric != "euclidean":
                logger.info("Falling back to euclidean distance")
                return self.sort_index(df, metric="euclidean", method=method)
            else:
                logger.info("Error encountered while sorting table:")
                self.log_df(df)
                raise e
        except Exception as e:
            logger.info("Error encountered while sorting table:")
            self.log_df(df)
            raise e

    def plot(self, meta_cname, output_folder):

        output_fp = f"{output_folder}/{meta_cname}"

        # Sort the bins and samples
        bin_order = self.sort_index(
            self.adata.to_df().T,
            metric="cosine",
            method="average"
        )

        sample_order = self.sort_index(
            self.adata.to_df(),
            metric="euclidean",
            method="ward"
        )

        genome_order = self.sort_index(
            self.adata.uns["group_profile"].T,
            metric="cosine",
            method="average"
        )

        # Relative size of rows and columns
        heatmap_size = 3
        column_widths = np.array([heatmap_size, 1, heatmap_size, 1, 1])
        column_widths = list(column_widths / column_widths.sum())
        row_heights = np.array([0.5, heatmap_size, 0.5])
        row_heights = list(row_heights / row_heights.sum())
        horizontal_spacing = 0.05

        cols = 5
        rows = 3
        fig = make_subplots(
            rows=rows,
            cols=cols,
            shared_xaxes=True,
            shared_yaxes=True,
            start_cell="bottom-left",
            column_widths=column_widths,
            row_heights=row_heights,
            horizontal_spacing=horizontal_spacing
        )
        # Bins across Genomes
        self.heatmap(
            (
                self
                .adata
                .uns["group_profile"]
                .reindex(
                    columns=genome_order,
                    index=bin_order
                )
            ),
            fig,
            value_label="Gene Bin Presence in Genome",
            row=2,
            col=1,
            showscale=True,
            coloraxis="coloraxis2",
            colorbar_x=(horizontal_spacing / 4),
            colorbar_xanchor="left",
            colorbar_orientation="h",
            colorbar_y=0.85,
            colorbar_len=column_widths[0] - (horizontal_spacing * 2),
            colorbar_title_side="top"
        )

        # # Bins Information

        # Number of genes per bin
        self.bar(
            pd.DataFrame(
                {
                    "Genes per Bin (#)": self.adata.uns["bin_size"]
                }
            ).reindex(
                index=bin_order
            ),
            "Genes per Bin (#)",
            fig,
            row=2,
            col=2,
            orient="h",
            log=True
        )

        # Bins across samples
        self.heatmap(
            self.log_scale((
                self
                .adata
                .to_df("prop")
                .reindex(
                    index=sample_order,
                    columns=bin_order
                )
                .T
            )),
            fig,
            value_label="Proportion of Reads (log10)",
            row=2,
            col=3,
            coloraxis="coloraxis3",
            showscale=True,
            colorbar_x=sum(column_widths[:2]) + (horizontal_spacing / 2),
            colorbar_xanchor="left",
            colorbar_orientation="h",
            colorbar_y=0.85,
            colorbar_len=column_widths[2] - (horizontal_spacing * 1.5),
            colorbar_title_side="top",
            colorbar_title_text="Gene Bin Abundance<br>Proportion of Reads (log10)",
            colorbar_xpad=0,
        )

        # Metadata data on the bins heatmap
        self.heatmap(
            self.adata.obs.reindex(
                index=sample_order,
                columns=[meta_cname]
            ).T,
            fig,
            value_label=meta_cname,
            col=3,
            row=1,
            colorscale=["blue", "orange"],
            coloraxis="coloraxis",
            showscale=False
        )

        # Bars showing the -log10(qvalue) for each gene bin
        # (unless there are no useful q-values, in which case fall back to p-values)
        if self.adata.var["neg_log10_qvalue"].max() < 0.1:
            kw = "neg_log10_pvalue"
            label = "p-value (-log10)"
        else:
            kw = "neg_log10_qvalue"
            label = "FDR-adjusted<br>p-value (-log10)"
        self.bar(
            (
                self
                .adata
                .var
                .reindex(
                    columns=[kw, "p_value", "estimate"]
                )
                .rename(
                    columns={kw: label}
                )
            ),
            label,
            fig,
            orient="h",
            col=4,
            row=2
        )

        # Boxplot comparing gene bin abundances by sample group
        self.box(
            self.adata.to_df("prop").reindex(index=sample_order),
            self.adata.obs[meta_cname].reindex(index=sample_order),
            meta_cname,
            "Gene Bin Abundance",
            fig,
            col=5,
            row=2,
            log=True,
            orient="h"
        )

        logger.info(fig.layout)
        title_text = f"Association with {meta_cname}"
        fig.update_layout(
            title_text=title_text,
            title_xanchor="center",
            title_x=0.5,
            plot_bgcolor='rgba(255, 255, 255, 1.0)',
            paper_bgcolor='rgba(255, 255, 255, 1.0)',
            legend=dict(
                yanchor="top",
                xanchor="left",
                x=1.0,
                y=0.3
            ),
            **{
                f"xaxis{i + cols}": dict(showticklabels=True)
                for i in [1, 2, 4, 5]
            },
            **{
                f"yaxis{2 + (cols * i)}": dict(showticklabels=True)
                for i in [1, 2, 4, 5, 6]
            }
        )
        fig.write_html(f"{output_fp}.html")
        for ext in ['png', 'pdf']:
            fig.write_image(
                f"{output_fp}.{ext}",
                width=2880,
                height=1800
            )

    @staticmethod
    def heatmap(
        df: pd.DataFrame,
        fig,
        value_label: str,
        row=None,
        col=None,
        coloraxis="coloraxis",
        colorscale=px.colors.sequential.Blues,
        showlegend=False,
        showscale=False,
        **kwargs
    ):

        assert df.dropna().shape[0] == df.shape[0], df

        if "colorbar_title_text" not in kwargs:
            kwargs["colorbar_title_text"] = value_label

        fig.layout[coloraxis] = dict(
            colorscale=colorscale,
            showscale=showscale,
            **{
                kw: val
                for kw, val in kwargs.items()
                if kw.startswith("colorbar_")
            }
        )

        fig.add_trace(
            go.Heatmap(
                z=df.values,
                x=df.columns.values,
                y=df.index.values,
                hovertext=df.applymap(
                    lambda v: f"{value_label}: {v}"
                ),
                coloraxis=coloraxis,
                showlegend=showlegend,
                showscale=showscale
            ),
            row=row,
            col=col
        )

    def bar(
        self,
        data: pd.DataFrame,
        cname: str,
        fig,
        row,
        col,
        orient="v",
        log=False
    ):
        assert orient in ["v", "h"]
        kwargs = (
            dict(
                y=data[cname].values,
                x=data[cname].index.values
            )
            if orient == "v" else
            dict(
                x=data[cname].values,
                y=data[cname].index.values,
                orientation='h'
            )
        )
        fig.add_trace(
            go.Bar(
                showlegend=False,
                hovertext=data.apply(
                    lambda r: "<br>".join([
                        f"{kw}: {val}"
                        for kw, val in r.items()
                    ]),
                    axis=1
                ),
                marker=dict(color="blue"),
                **kwargs
            ),
            row=row,
            col=col
        )
        self.label_axis(fig, row, col, orient, cname)
        if log:
            self.logscale_axis(fig, row, col, orient)

    def label_axis(self, fig, row, col, orient, title):
        self.update_axis(
            fig, row, col, orient, title=title
        )

    def logscale_axis(self, fig, row, col, orient):
        self.update_axis(
            fig, row, col, orient, type="log"
        )

    @staticmethod
    def update_axis(fig, row, col, orient, **kwargs):
        if orient == "v":
            fig.for_each_yaxis(
                lambda axis: axis.update(**kwargs),
                row=row,
                col=col
            )
        else:
            fig.for_each_xaxis(
                lambda axis: axis.update(**kwargs),
                row=row,
                col=col
            )

    def box(
        self,
        df: pd.DataFrame,
        huevals: pd.Series,
        huelabel: str,
        label: str,
        fig,
        row,
        col,
        orient="v",
        log=False,
        **box_kwargs
    ):
        assert orient in ["v", "h"]
        colors = ["blue", "orange", "red", "green", "black"]
        for ix, (hueval, hue_df) in enumerate(df.groupby(huevals)):
            trace_df = (
                hue_df
                .reset_index()
                .melt(
                    id_vars=[(
                        "index"
                        if hue_df.index.name is None
                        else hue_df.index.name
                    )],
                    var_name="variable"
                )
            )

            kwargs = (
                dict(
                    y=trace_df["value"],
                    x=trace_df["variable"]
                )
                if orient == "v" else
                dict(
                    orientation="h",
                    x=trace_df["value"],
                    y=trace_df["variable"]
                )
            )

            fig.add_trace(
                go.Box(
                    name=f"{huelabel}: {hueval}",
                    legendgroup=f"{huelabel}: {hueval}",
                    marker=dict(color=colors[ix % len(colors)]),
                    **kwargs,
                    **box_kwargs
                ),
                row=row,
                col=col
            )
        label = label.replace(" ", "<br>")
        self.label_axis(fig, row, col, orient, label)
        if log:
            self.logscale_axis(fig, row, col, orient)


def plot_metagenomes(
    param,
    stats,
    h5ad,
    output_folder
):
    mdata = Metagenome(stats, h5ad)
    mdata.plot(param, output_folder)


@click.command
@click.option('--param', type=str)
@click.option('--stats', type=click.Path(exists=True))
@click.option('--h5ad', type=click.Path(exists=True))
@click.option('--output_folder', type=click.Path())
def main(
    param,
    stats,
    h5ad,
    output_folder
):

    logger.info(f"param: {param}")
    logger.info(f"stats: {stats}")
    logger.info(f"counts: {h5ad}")
    logger.info(f"output_folder: {output_folder}")

    plot_metagenomes(
        param,
        stats,
        h5ad,
        output_folder
    )


if __name__ == "__main__":
    main()
