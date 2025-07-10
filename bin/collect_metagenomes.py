#!/usr/bin/env python3
from collections.abc import Mapping
import json
import os
import anndata as ad
import click
import logging
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_score
from mudata import MuData

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
                unaligned_reads=lambda d: d['tot_reads'] - d['aligned_reads'],
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

        # Remove any columns from obs which are all NaN
        obs = obs.loc[:, obs.notna().any()]

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
        min_n_genes,
        incl_unaligned
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
        self.calc_rel_abund(incl_unaligned != "false")

        return self

    def _read_csv(
        self,
        fp,
        index_col,
        cnames=[],
        fillna=False,
        **kwargs
    ) -> pd.DataFrame:

        logger.info(f"Reading in {fp}")
        df: pd.DataFrame = pd.read_csv(fp, **kwargs)

        for kw in [index_col] + cnames:
            msg = f"Expected to find {kw} in {fp}"
            assert kw in df.columns.values, msg

        # Drop any rows where the index_column is null
        df = df.dropna(subset=[index_col])

        # The index must be a string
        df = df.assign(**{index_col: df[index_col].astype(str)})
        df = df.set_index(index_col)

        # If fillna is True, fill NaN values with the provided value
        if fillna:
            df = df.fillna(fillna)

        return df

    def _read_csv_columns(
        self,
        fp=None,
        modality=None,
        index_col=None,
        attr=None,
        cnames=[],
        fillna=False,
        **kwargs
    ):
        """Add columns from CSV to <attr>"""

        for cname, cvals in self._read_csv(
            fp,
            index_col,
            cnames,
            fillna=fillna,
            **kwargs
        ).items():
            if cvals.dtype == "object":
                logger.info(f"Converting {cname} to a string")
                cvals = cvals.astype(str)
            logger.info(f"Adding {cname} to {attr} (type: {cvals.dtype})")
            getattr(
                (
                    self.data
                    if modality is None else
                    self.data[modality]
                ),
                attr
            )[(
                "_specimen"
                if cname == "specimen" and attr == "obs"
                else cname
            )] = cvals

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

        # Read in the table of genome groups
        genome_groups = self._read_csv(genome_groups, "genome")

        # Drop any columns which are entirely NaN
        genome_groups = genome_groups.loc[:, genome_groups.notna().any()]

        # Fill any remaining NaNs with "None"
        genome_groups = genome_groups.fillna("None")

        # For each column, if any item in the column is a string, convert
        # the column to a string
        for col in genome_groups.columns:
            if genome_groups[col].dtype == "object":
                genome_groups[col] = genome_groups[col].astype(str)

        # Save the list of which genomes belong to each group
        self.data.uns["genome_groups"] = genome_groups

        # Save the total number of genomes in each group
        self.data.uns["genome_group_size"] = (
            genome_groups
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
            cnames=[],
            fillna="None",
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
        logger.info(f"{filt_n:,} / {tot_n:,} samples filtered out: {query}")

        if "filtered_out" not in self.data.uns:
            self.data.uns["filtered_out"] = pd.Series(
                False,
                index=self.data.obs_names
            )

        self.data.uns["filtered_out"] = (
            self.data.uns["filtered_out"] |
            filtered_obs
        )
        filt_n = self.data.uns["filtered_out"].sum()
        logger.info(f"Total samples filtered out: {filt_n:,} / {tot_n:,}")

    @property
    def filtered_samples(self):
        """List of samples which pass the filters."""
        return self.data.uns["filtered_out"].index.values[
            ~self.data.uns["filtered_out"]
        ]

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
                .reindex(index=self.filtered_samples)
            )
        )
        self.data.mod["bins"].var["n_genes"] = pd.Series(self.data.uns["bin_size"])
        self.data.update()
        self.log_content()

    def calc_bin_silhouette_score(self, metric="euclidean"):
        """Calculate the silhouette score for each bin."""

        # If there are more than 10,000 genes, don't calculate
        if self.data.mod["genes"].to_df().shape[1] > 10000:
            logger.info("Too many genes to calculate silhouette score")
            return

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
        group_profile = (group_profile / group_profile.sum()).fillna(0)

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
            # Run NNLS, skipping if there are any errors
            logger.info(f"Attempting NNLS on sample {sample}")
            try:
                x, rnorm = nnls(
                    group_profile,
                    bin_abund
                )
            except (ValueError, RuntimeError, np.linalg.LinAlgError):
                logger.info(f"Sample failed to run through NNLS - skipping ({sample})")
                continue

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
                    index=self.filtered_samples,
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
                    index=self.filtered_samples
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

    def calc_rel_abund(self, incl_unaligned: bool, mods=None):
        """
        Normalize the abundance to the number of reads aligned

        For each modality, a 'prop' layer will be created
        """
        if mods is None:
            mods = ["bins", "genomes"]

        for mod in mods:
            # Get the values for this modality
            abund: pd.DataFrame = self.data.mod[mod].to_df()

            assert "aligned_reads" in self.data.mod["genes"].obs
            logger.info(f"Normalizing {mod} to sequencing depth")
            rel_abund = (
                abund.T /
                self.data.mod["genes"].obs[
                    "tot_reads" if incl_unaligned else "aligned_reads"
                ].reindex(
                    index=abund.index.values
                )
            ).T
            self.data.mod[mod].layers["prop"] = rel_abund

    def to_h5ad(self, output_folder):

        os.makedirs(output_folder, exist_ok=True)
        if "filtered_out" in self.data.uns:
            del self.data.uns["filtered_out"]
        for mod in self.data.mod:
            adata: ad.AnnData = self.data.mod[mod]
            adata.obs = (
                self.data
                .obs
                .reindex(index=adata.obs.index)
            )
            # Remove any columns from obs which are all NaN
            adata.obs = adata.obs.loc[:, adata.obs.notna().any()]
            # Fill any missing values in obs with "None"
            adata.obs = adata.obs.fillna("None")

            for kw, val in self.data.uns.items():
                # Omit the group profile from the uns
                if kw != "group_profile":
                    adata.uns[kw] = val
            for line in str(adata).split("\n"):
                logger.info(line)
            logger.info(f"Writing to: {output_folder}/metagenome.{mod}.h5ad")
            adata.write_h5ad(f"{output_folder}/metagenome.{mod}.h5ad", compression="gzip")

    def to_csv(self, output_folder):

        os.makedirs(output_folder, exist_ok=True)
        if "filtered_out" in self.data.uns:
            del self.data.uns["filtered_out"]
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

    @staticmethod
    def log_scale(df: pd.DataFrame):

        lowest = df.apply(lambda c: c[c > 0].min()).min()
        return df.clip(lower=lowest).apply(np.log10)

    def sort_index(self, df: pd.DataFrame, metric="cosine", method="average"):

        if df.shape[0] < 3:
            logger.info(f"No need to sort - {df.shape[0]} rows")
            return df.index.values
        logger.info(f"Sorting {df.shape[0]:,} rows")

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
        except Exception as e:
            logger.info("Error encountered while sorting table:")
            self.log_df(df)
            raise e


def bin_metagenomes(
    read_alignments,
    gene_bins,
    genome_groups,
    group_profile,
    metadata,
    incl_unaligned,
    min_n_reads,
    min_n_genes
):
    mdata = Metagenome.from_read_alignments(
        read_alignments,
        gene_bins,
        genome_groups,
        group_profile,
        min_n_reads,
        min_n_genes,
        incl_unaligned
    )
    mdata.add_metadata(metadata)
    return mdata


@click.command
@click.option('--read_alignments', type=click.Path(exists=True))
@click.option('--gene_bins', type=click.Path(exists=True))
@click.option('--genome_groups', type=click.Path(exists=True))
@click.option('--group_profile', type=click.Path(exists=True))
@click.option('--metadata', type=click.Path(exists=True))
@click.option('--incl_unaligned', type=str)
@click.option('--min_n_reads', type=int)
@click.option('--min_n_genes', type=int)
@click.option('--output_folder', type=click.Path())
def main(
    read_alignments,
    gene_bins,
    genome_groups,
    group_profile,
    metadata,
    incl_unaligned,
    min_n_reads,
    min_n_genes,
    output_folder
):

    logger.info(f"read_alignments: {read_alignments}")
    logger.info(f"gene_bins: {gene_bins}")
    logger.info(f"genome_groups: {genome_groups}")
    logger.info(f"group_profile: {group_profile}")
    logger.info(f"metadata: {metadata}")
    logger.info(f"incl_unaligned: {incl_unaligned}")
    logger.info(f"min_n_reads: {min_n_reads}")
    logger.info(f"min_n_genes: {min_n_genes}")

    mdata = bin_metagenomes(
        read_alignments,
        gene_bins,
        genome_groups,
        group_profile,
        metadata,
        incl_unaligned,
        min_n_reads,
        min_n_genes
    )

    mdata.to_csv(f"{output_folder}/csv")
    mdata.to_h5ad(f"{output_folder}/h5ad")
    mdata.data.write_h5mu(f"{output_folder}/metagenome.h5mu")


if __name__ == "__main__":
    main()
