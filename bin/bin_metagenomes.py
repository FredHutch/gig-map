#!/usr/bin/env python3
from collections.abc import Mapping
import json
import anndata as ad
import click
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy import stats
from scipy.optimize import nnls
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests
from muon import MuData

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [bin_metagenomes.py] %(message)s'
)
logger = logging.getLogger('bin_metagenomes.py')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Write to file
fileHandler = logging.FileHandler("bin_metagenomes.log")
fileHandler.setFormatter(logFormatter)
fileHandler.setLevel(logging.INFO)
logger.addHandler(fileHandler)


class Metagenome:

    mdata: MuData
    fdr_method = "fdr_bh"

    """
    The mdata attribute is a MuData object with modalities
    that will include: genes, gene bins, and genomes

    Observations are specimens (i.e. samples)
    MuData object:
        obs:	'group'
        var:	'length'
        uns:	'bin_size', 'genome_group_size', 'group_profile'
        3 modalities
            genes:
                obs:	'tot_reads', 'aligned_reads', 'detected_genes'
                var:	'length', 'combined_name', 'bin', 'n_genomes'
                layers:	'std', 'depth', 'coverage'
                varp:	'pdist_euclidean'
            bins:
                var:	'n_genes', 'silhouette_score'
                varm:	'~ category - Relative to Bin 1', '~ category - Relative to Total Reads', '~ category - Relative to Aligned Reads'
                layers:	'Relative to Bin 1', 'Relative to Total Reads', 'Relative to Aligned Reads'
            genomes:
                obs:	'nnls_residual'
                var:	'n_genomes'
                varm:	'~ category - Relative to Total Reads', '~ category - Relative to Aligned Reads'
                layers:	'Relative to Total Reads', 'Relative to Aligned Reads'
    """

    def __init__(self, **data):
        self.data = MuData(data)

    @classmethod
    def parse_read_alignments(cls, read_alignments_fp) -> 'Metagenome':

        logger.info(f"Reading in {read_alignments_fp}")
        read_alignments = pd.read_csv(read_alignments_fp)

        for kw in [
            "std",
            "depth",
            "length",
            "coverage",
            "nreads",
            "id",
            "specimen",
            "tot_reads"
        ]:
            msg = f"Expected to find {kw} in {read_alignments_fp}"
            assert kw in read_alignments.columns.values, msg

        obs_ix = read_alignments["specimen"].unique()
        var_ix = read_alignments["id"].unique()

        layers = {
            kw: (
                read_alignments
                .pivot(
                    index="specimen",
                    columns="id",
                    values=kw
                )
                .fillna(0)
                .reindex(
                    index=obs_ix,
                    columns=var_ix
                )
            )
            for kw in ["std", "depth", "coverage", "nreads"]
        }
        X = layers.pop("nreads")

        obs = (
            read_alignments
            .reindex(
                columns=["specimen", "tot_reads"]
            )
            .drop_duplicates()
            .set_index("specimen")
            .reindex(
                index=obs_ix
            )
            .assign(
                aligned_reads=X.sum(axis=1),
                detected_genes=(X > 0).sum(axis=1)
            )
        )

        var = (
            read_alignments
            .reindex(
                columns=["id", "length"]
            )
            .drop_duplicates()
            .set_index("id")
            .reindex(
                index=var_ix
            )
        )

        self = cls(
            genes=ad.AnnData(
                X=X,
                obs=obs,
                var=var,
                layers=layers
            )
        )
        self.data.update()
        return self

    @classmethod
    def from_read_alignments(
        cls,
        read_alignments_fp,
        gene_bins,
        genome_groups,
        group_profile,
        min_n_reads,
        min_n_genes
    ):
        self = cls.parse_read_alignments(read_alignments_fp)
        self.log_content()
        self.add_gene_bins(gene_bins)
        self.add_genome_groups(genome_groups)
        self.add_group_profile(group_profile)
        self.filter_by_aligned_reads(min_n_reads)
        self.filter_by_aligned_genes(min_n_genes)
        self.calc_bin_abund()
        self.calc_bin_silhouette_score()
        self.est_genome_abund()
        self.calc_rel_abund()

        return self

    def _read_csv(
        self,
        fp,
        index_col,
        cnames=[],
        **kwargs
    ) -> pd.DataFrame:

        logger.info(f"Reading in {fp}")
        df: pd.DataFrame = pd.read_csv(fp, **kwargs)

        for kw in [index_col] + cnames:
            msg = f"Expected to find {kw} in {fp}"
            assert kw in df.columns.values, msg

        return df.set_index(index_col)

    def _read_csv_columns(
        self,
        fp=None,
        modality=None,
        index_col=None,
        attr=None,
        cnames=[],
        **kwargs
    ):
        """Add columns from CSV to <attr>"""

        for cname, cvals in self._read_csv(
            fp,
            index_col,
            cnames,
            **kwargs
        ).items():
            logger.info(f"Adding {cname} to {attr}")
            getattr(
                (
                    self.data
                    if modality is None else
                    self.data[modality]
                ),
                attr
            )[cname] = cvals

        self.log_content()

    def add_gene_bins(self, gene_bins):
        logger.info("Adding gene bins")

        # Add the bin annotation to var
        self._read_csv_columns(
            fp=gene_bins,
            modality="genes",
            index_col="gene_id",
            attr="var",
            cnames=["bin"]
        )
        # Save the total number of genes in each bin
        # (even if they aren't found in the alignments)
        self.data.uns["bin_size"] = (
            self._read_csv(
                gene_bins,
                "gene_id"
            )
            ["bin"]
            .value_counts()
            .to_dict()
        )

    def add_genome_groups(self, genome_groups):
        logger.info("Adding genome groups")

        # Save the total number of genomes in each group
        self.data.uns["genome_group_size"] = (
            self._read_csv(
                genome_groups,
                "genome"
            )
            ["group"]
            .value_counts()
            .to_dict()
        )

    def add_group_profile(self, group_profile):
        logger.info("Adding group profile")

        # Save the proportion of genes from each bin
        # which are found in each group
        self.data.uns["group_profile"] = (
            self._read_csv(
                group_profile,
                "group"
            )
            .T
            .reindex(
                index=pd.Series(self.data.uns["bin_size"]).index
            )
        )
        self.log_content()

    def add_metadata(self, metadata):
        logger.info("Adding metadata")
        self._read_csv_columns(
            fp=metadata,
            index_col="sample",
            attr="obs",
            cnames=[]
        )
        self.data.update()
        self.log_content()

    def filter_by_aligned_reads(self, min_n_reads) -> 'Metagenome':
        """Filter specimens based on a minimum aligned read depth."""
        logger.info(f"Filtering by >={min_n_reads} aligned reads")

        self._filter_self(
            (self.data.mod["genes"].obs["aligned_reads"] < min_n_reads),
            f"< {min_n_reads:,} reads aligned"
        )

    def filter_by_aligned_genes(self, min_n_genes) -> 'Metagenome':
        """Filter specimens based on a minimum number of genes detected."""
        logger.info(f"Filtering by >={min_n_genes} aligned genes")
        self._filter_self(
            (self.data.mod["genes"].obs["detected_genes"] < min_n_genes),
            f"< {min_n_genes:,} genes detected"
        )

    def _filter_self(
        self,
        filtered_obs: pd.Series,
        query: str
    ) -> 'Metagenome':

        filt_n = filtered_obs.sum()
        tot_n = filtered_obs.shape[0]
        logger.info(f"{filt_n:,} / {tot_n:,} samples have {query}")

        if filtered_obs.any():
            self.data = self.data[~filtered_obs]
        self.log_content()

    def log_content(self):
        for line in str(self.data).split("\n"):
            logger.info(line)

    @staticmethod
    def log_df(df: pd.DataFrame, **kwargs):
        for line in df.to_csv(**kwargs).split("\n"):
            logger.info(line)

    def calc_bin_abund(self):
        """
        Calculate the abundance (number of reads) from each gene bin
        across every sample.
        """
        logger.info("""Calculating the sequencing depth per bin""")

        self.data.mod["bins"] = ad.AnnData(
            (
                self.data.mod["genes"]
                .to_df()
                .T
                .groupby(
                    self.data.mod["genes"].var["bin"]
                )
                .sum()
                .T
            )
        )
        self.data.mod["bins"].var["n_genes"] = pd.Series(self.data.uns["bin_size"])
        self.data.update()
        self.log_content()

    def calc_bin_silhouette_score(self, metric="euclidean"):
        """Calculate the silhouette score for each bin."""
        logger.info("Calculating silhouette score for bins")

        # Calculate pairwise distances
        self.calc_pdist_var("genes", metric)

        scores = {}
        for bin in self.data.mod["genes"].var["bin"].dropna().unique():

            labels = (self.data.mod["genes"].var["bin"] == bin).apply(int)
            scores[bin] = silhouette_score(
                self.data.mod["genes"].varp[f"pdist_{metric}"],
                labels,
                metric="precomputed"
            )

        self.data.mod["bins"].var["silhouette_score"] = (
            pd.Series(scores)
            .reindex(
                index=self.data.mod["bins"].var_names
            )
        )
        self.log_content()

    def calc_pdist_var(self, mod, metric="euclidean"):
        logger.info(f"Calculating {metric} distances for {mod} variables")

        # Calculate the distance matrix for all genes
        kw = f"pdist_{metric}"
        self.data.mod[mod].varp[kw] = distance.squareform(
            distance.pdist(
                self.data.mod[mod].to_df().T.values,
                metric=metric
            )
        )
        self.log_content()

    def est_genome_abund(self):
        """
        Estimate the relative abundance of each genome group
        using a mixture modeling approach.
        """
        logger.info("Estimating genome abundances")

        # Keep track of the NNLS residual for each sample
        nnls_residual = {}

        # Build a dict with the estimated abundances for each genome group
        genome_abund = {}

        # Scale the group profile by the number of genes in each bin
        group_profile: pd.DataFrame = (
            self.data.uns["group_profile"].T *
            pd.Series(self.data.uns["bin_size"])
        ).T

        # Translate the scaled gene abundance into a proportion, so
        # that the fitted data will ideally have the same sum as
        # the number of reads sequenced for that sample
        group_profile = group_profile / group_profile.sum()

        # Make a table with the number of reads per bin
        reads_per_bin = (
            self.data
            .mod["bins"]
            .to_df()
            .reindex(
                columns=group_profile.index.values
            )
            .fillna(0)
        )

        for sample, bin_abund in reads_per_bin.iterrows():
            # Run NNLS
            x, rnorm = nnls(
                group_profile,
                bin_abund
            )

            # Scale the 2-norm of the residual by the total
            # number of reads in the right-hand side vector
            nnls_residual[sample] = rnorm / bin_abund.sum()

            # Add a vector with the fitted genome abundances
            genome_abund[sample] = (
                pd.Series(
                    x,
                    index=group_profile.columns.values
                )
            )

        # Make a DataFrame with the fitted read numbers per genome
        genome_abund = pd.DataFrame(genome_abund)

        self.data.mod["genomes"] = ad.AnnData(
            X=(
                genome_abund
                .T
                .reindex(
                    index=self.data.obs_names,
                    columns=group_profile.columns.values
                )
            ),
            obs=(
                pd.DataFrame(
                    dict(
                        nnls_residual=nnls_residual
                    )
                )
                .reindex(
                    index=self.data.obs_names
                )
            ),
            var=(
                pd.DataFrame(
                    dict(
                        n_genomes=self.data.uns["genome_group_size"]
                    )
                )
                .reindex(
                    index=group_profile.columns.values
                )
            )
        )
        self.data.update()
        self.log_content()

    def calc_rel_abund(self, mods=None, prefix="Relative to "):
        """
        Compute the relative abundance in different ways.
        1. Relative to any feature which is ubiquitously present
        2. Relative to the total number of reads sequenced
        3. Relative to the number of reads aligned

        For each modality, a layer will be created with the given prefix
        """
        if mods is None:
            mods = ["bins", "genomes"]

        for mod in mods:
            # Get the values for this modality
            abund: pd.DataFrame = self.data.mod[mod].to_df()

            # Identify any variables which are found in all samples
            comp_vars = [var for var, vals in abund.items() if (vals > 0).all()]
            logger.info(f"There are {len(comp_vars):,} variables {mod} present in all samples")
            for var in comp_vars:
                kw = f"{prefix}{var}"
                logger.info(f"Computing {mod} {kw}")
                # Calculate the ratio to the comparison variable
                ratios = (abund.T / abund[var]).T
                # Save the table
                self.data.mod[mod].layers[kw] = ratios

            # Run the analysis in comparison to the number of aligned reads
            for kw, label in [
                ("tot_reads", "Total Reads"),
                ("aligned_reads", "Aligned Reads")
            ]:
                assert kw in self.data.mod["genes"].obs
                new_kw = f"{prefix}{label}"
                logger.info(f"Computing {mod} {new_kw}")
                rel_abund = (abund.T / self.data.mod["genes"].obs[kw]).T
                self.data.mod[mod].layers[new_kw] = rel_abund

    def compare_groups(self, category, mods=None):
        assert category in self.data.obs.columns.values

        # Get the metadata values
        meta: pd.Series = self.data.obs[category].dropna()
        logger.info(f"Samples with {category} defined: {meta.shape[0]:,}")

        # We can only support the comparison of two groups (at this point)
        logger.info(f"Unique values: {', '.join([str(v) for v in meta.unique()])}")
        assert meta.unique().shape[0] == 2, f"{category} must contain 2 groups"

        # Convert meta to 1/0
        meta = (meta == max(meta.unique())).apply(int)

        if mods is None:
            mods = ["bins", "genomes"]

        for mod in mods:
            self.compare_groups_mod(meta, mod)
            self.log_content()

    @staticmethod
    def pack_meta_layer(meta_name, layer):
        return f"~ {meta_name} - {layer}"

    @staticmethod
    def unpack_meta_layer(varm):
        if varm.startswith("~ "):
            return varm[2:].split(" - ", 1)
        else:
            return None, None

    def compare_groups_mod(self, meta: pd.Series, mod: str):
        logger.info(f"Comparing samples on the basis of {mod}")

        # Use each of the different relative abundance tables
        # which have been computed for this modality
        for layer in self.data.mod[mod].layers:
            label = self.pack_meta_layer(meta.name, layer)
            logger.info(f"Comparing samples {label}")
            abund: pd.DataFrame = self.data.mod[mod].to_df(layer)

            self.data.mod[mod].varm[label] = (
                self.compare_groups_single(meta, abund)
            )

    def compare_groups_single(
        self,
        meta: pd.Series,
        abund: pd.DataFrame
    ):
        # At this point meta should just be 0 or 1
        assert meta.isin([0, 1]).all()
        comparison_obs = meta.index.values[meta == 0]
        control_obs = meta.index.values[meta == 1]

        res = (
            pd.DataFrame([
                self.mannwhitneyu(var, vals, control_obs, comparison_obs)
                for var, vals in abund.items()
            ])
            .set_index("index")
            .reindex(index=abund.columns.values)
        )

        res = res.assign(
            qvalue=multipletests(res["pvalue"].values, method=self.fdr_method)[1],
            neg_log10_pvalue=lambda d: d['pvalue'].apply(np.log10) * -1,
            neg_log10_qvalue=lambda d: d['qvalue'].apply(np.log10) * -1
        )
        return res

    @staticmethod
    def mannwhitneyu(var: str, vals: pd.Series, control_obs, comparison_obs):
        res = stats.mannwhitneyu(
            vals.reindex(index=control_obs).values,
            vals.reindex(index=comparison_obs).values
        )
        return dict(
            pvalue=res.pvalue,
            statistic=res.statistic,
            index=var
        )

    def to_h5ad(self, output_folder):

        for mod in self.data.mod:
            self.data.mod[mod].write_h5ad(f"{output_folder}/metagenome.{mod}.h5ad")

    def to_csv(self, output_folder):

        self.to_csv_obj(self.data, f"{output_folder}/metagenome")

    def to_csv_obj(self, obj, path):

        if isinstance(obj, MuData):

            if obj.obs.shape[1] > 0:
                obj.obs.to_csv(path + ".obs.csv.gz")
            for kw, val in obj.uns.items():
                self.to_csv_obj(val, path + "." + kw)
            for kw, val in obj.mod.items():
                self.to_csv_obj(val, path + "." + kw)

        elif isinstance(obj, (dict, list, str, int, float)):
            self.write_json(obj, path + ".json")

        elif isinstance(obj, pd.DataFrame):
            obj.to_csv(path + ".csv.gz")

        elif isinstance(obj, ad.AnnData):
            for kw in dir(obj):
                if kw.startswith(('obs', 'var')):
                    self.to_csv_obj(
                        getattr(obj, kw),
                        path + "." + kw
                    )
            for layer in obj.layers:
                self.to_csv_obj(
                    obj.to_df(layer),
                    path + "." + layer
                )
            self.to_csv_obj(obj.to_df(), path + '.X')

        elif isinstance(obj, Mapping):
            for kw, val in obj.items():
                self.to_csv_obj(
                    val,
                    path + "." + kw
                )

    @staticmethod
    def write_json(val, path):
        if isinstance(val, pd.Series):
            val = val.to_dict()
        with open(path, "w") as handle:
            json.dump(val, handle, indent=4)

    def plot(self, output_folder):

        # The central heatmap is the comparison of bins and genomes
        # Bins extend vertically from that heatmap,
        # and genomes extend horizontally.

        for bins_varm in self.data.mod["bins"].varm:
            bins_meta, bins_layer = self.unpack_meta_layer(bins_varm)

            for genomes_varm in self.data.mod["genomes"].varm:
                genomes_meta, genome_layer = self.unpack_meta_layer(genomes_varm)

                if bins_meta == genomes_meta:
                    output_fn = ' - '.join([bins_meta, "bins", bins_layer, "genomes", genome_layer])
                    output_fp = f"{output_folder}/{output_fn}"
                    self.write_image(
                        bins_meta,
                        bins_layer,
                        genome_layer,
                        output_fp,
                    )

    @staticmethod
    def sort_index(df):
        return df.index.values[
            hierarchy.leaves_list(
                hierarchy.linkage(
                    df.values,
                    metric="cosine",
                    method="average"
                )
            )
        ]

    def write_image(
        self,
        meta_cname,
        bins_layer,
        genome_layer,
        output_fp
    ):
        """"
        Figure layout:
        Genomes vs. Bins in the center, with bins extending
        vertically and genomes extending horizontally.

        |                         | (2, 7): Bins Boxplot    |
        |                         | (2, 6): Bins -log10(p)  |
        |                         | (2, 5): Bins Silhouette |
        | (1, 3): Sample Metadata | (2, 4): Bins Heatmap    | 
        |                         | (2, 3): # Genes per Bin |                      | (4, 3): NNLS Residual   |
        |                         | (2, 2): Central Heatmap | (3, 2): # of Genomes | (4, 2): Genomes Heatmap | (5, 2): Genomes -log10(p) | (6, 2): Genomes Boxplot |
        |                         |                         |                      | (4, 1): Sample Metadata |
        """

        # Relative size of rows and columns
        heatmap_size = 3
        column_widths = np.array([0.5, heatmap_size, 1, heatmap_size, 1, 1])
        column_widths = list(column_widths / column_widths.sum())
        row_heights = np.array([0.5, heatmap_size, 1, heatmap_size, 1, 1, 1])
        row_heights = list(row_heights / row_heights.sum())

        # Sort the bins and genomes
        bin_order = self.sort_index(self.data.uns["group_profile"])
        genome_order = self.sort_index(self.data.uns["group_profile"].T)

        cols = 6
        rows = 7
        fig = make_subplots(
            rows=rows,
            cols=cols,
            shared_xaxes=True,
            shared_yaxes=True,
            start_cell="bottom-left",
            column_widths=column_widths,
            row_heights=row_heights
        )
        # Bins across Genomes
        self.heatmap(
            (
                self.data
                .uns["group_profile"]
                .T
                .reindex(
                    index=genome_order,
                    columns=bin_order
                )
            ),
            fig,
            value_label="Gene Bin Presence in Genome",
            row=2,
            col=2,
            showscale=True,
            coloraxis="coloraxis2",
            colorbar_x=0.0,
            colorbar_xpad=0,
            colorbar_title_side="right",
            colorbar_xanchor="center",
            colorbar_yanchor="top",
            colorbar_len=0.25,
            colorbar_y=0.3
        )

        # Number of genomes per group
        self.bar(
            (
                self.data
                .mod["genomes"]
                .var["n_genomes"]
                .reindex(index=genome_order)
            ),
            "Number of Genomes per Group",
            fig,
            row=2,
            col=3,
            orient="h"
        )

        # Genomes across samples
        genomes_df: pd.DataFrame = (
            self.data
            .mod["genomes"]
            .to_df(genome_layer)
        )
        sample_order = self.sort_index(genomes_df)
        genomes_df = genomes_df.reindex(
            columns=genome_order,
            index=sample_order
        ).fillna(0)

        if genome_layer.startswith("Relative to "):
            comparitor = genome_layer[len("Relative to "):]
            if comparitor in genomes_df.columns.values:
                genomes_df = (
                    genomes_df
                    .drop(columns=[comparitor])
                    .reindex(columns=genome_order)
                )

        self.heatmap(
            genomes_df.T,
            fig,
            value_label=genome_layer,
            row=2,
            col=4,
            showscale=True,
            coloraxis="coloraxis3",
            colorbar_x=0.5,
            colorbar_xanchor="left",
            colorbar_xpad=0,
            colorbar_orientation="h",
            colorbar_title_side="top",
            colorbar_title_text=f"Genome Group Abundance<br>{genome_layer}",
            colorbar_y=0.45,
            colorbar_len=0.25
        )

        # Annotate metadata on the genomes heatmap
        metadata = (
            self.data
            .obs
            .reindex(
                columns=[meta_cname],
                index=sample_order
            )
        )

        self.heatmap(
            metadata.reindex(index=sample_order).T,
            fig,
            value_label=meta_cname,
            row=1,
            col=4,
            colorscale=["blue", "orange"],
            showscale=False
        )

        # NNLS Residual on genome abundance predictions from samples
        self.bar(
            (
                self.data
                .mod["genomes"]
                .obs["nnls_residual"]
                .reindex(index=sample_order)
            ),
            "NNLS Residual",
            fig,
            row=3,
            col=4
        )
        fig.for_each_yaxis(
            lambda axis: axis.update(matches=None, showticklabels=True),
            row=3,
            col=4
        )

        # Bars showing the -log10(qvalue) for each genome group
        self.bar(
            (
                self.data
                .mod["genomes"]
                .varm[f"~ {meta_cname} - {genome_layer}"]
                ["neg_log10_qvalue"]
            ),
            "FDR-adjusted p-value (-log10)",
            fig,
            row=2,
            col=5,
            orient="h"
        )

        # Boxplot comparing genome abundances by sample group
        self.box(
            genomes_df,
            self.data.obs[meta_cname].reindex(index=sample_order),
            meta_cname,
            "Genome Group Abundance",
            fig,
            row=2,
            col=6,
            orient="h",
            showlegend=False
        )

        # Bins Information (extending vertically from the central heatmap)

        # Number of genes per bin
        self.bar(
            pd.Series(
                self.data.uns["bin_size"]
            ).reindex(
                index=bin_order
            ),
            "Number of Genes per Bin",
            fig,
            row=3,
            col=2
        )

        # Bins across samples
        bins_df: pd.DataFrame = (
            self.data
            .mod["bins"]
            .to_df(bins_layer)
        )
        sample_order = self.sort_index(bins_df)
        bins_df = bins_df.reindex(
            columns=bin_order,
            index=sample_order
        ).fillna(0)
        if bins_layer.startswith("Relative to "):
            comparitor = bins_layer[len("Relative to "):]
            if comparitor in bins_df.columns.values:
                bins_df = (
                    bins_df
                    .drop(columns=[comparitor])
                    .reindex(columns=bin_order)
                )

        self.heatmap(
            bins_df,
            fig,
            value_label=bins_layer,
            row=4,
            col=2,
            coloraxis="coloraxis4",
            showscale=True,
            colorbar_x=0.35,
            colorbar_xanchor="left",
            colorbar_y=0.55,
            colorbar_len=0.25,
            colorbar_title_side="right",
            colorbar_title_text=f"Gene Bin Abundance<br>{bins_layer}",
            colorbar_xpad=0
        )

        # Metadata data on the bins heatmap
        self.heatmap(
            metadata.reindex(index=sample_order),
            fig,
            value_label=meta_cname,
            row=4,
            col=1,
            colorscale=["blue", "orange"],
            coloraxis="coloraxis",
            showscale=False
        )

        # Silhouette score on bin assignments
        self.bar(
            (
                self.data
                .mod["bins"]
                .var["silhouette_score"]
                .reindex(index=bin_order)
            ),
            "Silhouette Score",
            fig,
            row=5,
            col=2
        )

        # Bars showing the -log10(qvalue) for each gene bin
        self.bar(
            (
                self.data
                .mod["bins"]
                .varm[f"~ {meta_cname} - {bins_layer}"]
                ["neg_log10_qvalue"]
            ),
            "FDR-adjusted p-value (-log10)",
            fig,
            row=6,
            col=2
        )

        # Boxplot comparing genome abundances by sample group
        self.box(
            bins_df,
            self.data.obs[meta_cname].reindex(index=sample_order),
            meta_cname,
            "Gene Bin Abundance",
            fig,
            row=7,
            col=2
        )

        logger.info(fig.layout)
        title_text = " - ".join(
            [
                f"Association with {meta_cname}",
                f"Gene Bins: {bins_layer}",
                f"Genome Groups: {genome_layer}"
            ]
        )
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
                for i in [2, 3, 5, 6]
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

    def bar(self, vals: pd.Series, label: str, fig, row, col, orient="v"):
        assert orient in ["v", "h"]
        kwargs = (
            dict(
                y=vals.values,
                x=vals.index.values
            )
            if orient == "v" else
            dict(
                x=vals.values,
                y=vals.index.values,
                orientation='h'
            )
        )
        fig.add_trace(
            go.Bar(
                showlegend=False,
                hovertext=vals.apply(
                    lambda v: f"{label}: {v}"
                ),
                marker=dict(color="blue"),
                **kwargs
            ),
            row=row,
            col=col
        )
        label = label.replace(" ", "<br>")
        self.label_axis(fig, row, col, orient, label)

    @staticmethod
    def label_axis(fig, row, col, orient, label):
        if orient == "v":
            fig.for_each_yaxis(
                lambda axis: axis.update(title=label),
                row=row,
                col=col
            )
        else:
            fig.for_each_xaxis(
                lambda axis: axis.update(title=label),
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
        **box_kwargs
    ):
        assert orient in ["v", "h"]
        colors = ["blue", "orange", "red", "green", "black"]
        for ix, (hueval, hue_df) in enumerate(df.groupby(huevals)):
            trace_df = (
                hue_df
                .reset_index()
                .melt(
                    id_vars=[hue_df.index.name],
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


def bin_metagenomes(
    read_alignments,
    gene_bins,
    genome_groups,
    group_profile,
    metadata,
    category,
    min_n_reads,
    min_n_genes
):
    mdata = Metagenome.from_read_alignments(
        read_alignments,
        gene_bins,
        genome_groups,
        group_profile,
        min_n_reads,
        min_n_genes
    )
    mdata.add_metadata(metadata)
    mdata.compare_groups(category)
    return mdata


@click.command
@click.option('--read_alignments', type=click.Path(exists=True))
@click.option('--gene_bins', type=click.Path(exists=True))
@click.option('--genome_groups', type=click.Path(exists=True))
@click.option('--group_profile', type=click.Path(exists=True))
@click.option('--metadata', type=click.Path(exists=True))
@click.option('--category', type=str)
@click.option('--min_n_reads', type=int)
@click.option('--min_n_genes', type=int)
@click.option('--output_folder', type=click.Path())
def main(
    read_alignments,
    gene_bins,
    genome_groups,
    group_profile,
    metadata,
    category,
    min_n_reads,
    min_n_genes,
    output_folder
):

    logger.info(f"read_alignments: {read_alignments}")
    logger.info(f"gene_bins: {gene_bins}")
    logger.info(f"genome_groups: {genome_groups}")
    logger.info(f"group_profile: {group_profile}")
    logger.info(f"metadata: {metadata}")
    logger.info(f"category: {category}")
    logger.info(f"min_n_reads: {min_n_reads}")
    logger.info(f"min_n_genes: {min_n_genes}")

    mdata = bin_metagenomes(
        read_alignments,
        gene_bins,
        genome_groups,
        group_profile,
        metadata,
        category,
        min_n_reads,
        min_n_genes
    )

    # FIXME
    # mdata.to_csv(output_folder)
    # mdata.to_h5ad(output_folder)
    # mdata.data.write_h5mu(f"{output_folder}/metagenome.h5mu")
    mdata.plot(output_folder)


if __name__ == "__main__":
    main()
