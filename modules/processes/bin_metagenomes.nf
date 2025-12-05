process bin_summary {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path read_alignments
    path gene_bins
    path centroids_length

    output:
    path "bin_summary.csv.gz", emit: bin_summary

    """#!/bin/bash
set -e

bin_summary.py \
    --read-alignments "${read_alignments}" \
    --gene-bins "${gene_bins}" \
    --centroids-length "${centroids_length}"
    """
}

process wide_bin_abundance {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/bin_abundance", mode: 'copy', overwrite: true

    input:
    path "bin_summary.csv.gz"

    output:
    path "rpkm.csv.gz", emit: rpkm
    path "n_reads_aligned.csv.gz", emit: n_reads_aligned
    path "prop_reads_aligned.csv.gz", emit: prop_reads_aligned
    path "fragments_per_million.csv.gz", emit: fragments_per_million

    """#!/bin/bash
set -e

wide_bin_abundance.py \
    --min-n-reads ${params.min_n_reads} \
    --min-n-genes ${params.min_n_genes}
    """
}

process collect {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path read_alignments
    path gene_bins
    path genome_groups
    path group_profile
    path metadata

    output:
    path "*", emit: all
    path "csv/metagenome.obs.csv.gz", emit: metadata
    path "csv/metagenome.bins.X.csv.gz", emit: bin_counts
    path "h5ad/metagenome.bins.h5ad", emit: bins_h5ad
    path "h5ad/metagenome.genes.h5ad", emit: genes_h5ad
    path "h5ad/metagenome.genomes.h5ad", emit: genomes_h5ad


    """#!/bin/bash
set -e

MPLCONFIGDIR=\${TMPDIR:-/tmp/.config/matplotlib} \
collect_metagenomes.py \
    --read_alignments "${read_alignments}" \
    --gene_bins "${gene_bins}" \
    --genome_groups "${genome_groups}" \
    --group_profile "${group_profile}" \
    --metadata "${metadata}" \
    --incl_unaligned ${params.incl_unaligned} \
    --min_n_reads ${params.min_n_reads} \
    --min_n_genes ${params.min_n_genes} \
    --output_folder ./
    """
}

// Split up the regress results as inputs for plotting
process split {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
    path "regress.results.csv"

    output:
    path "*.results.csv"

    """#!/usr/bin/env python3
import pandas as pd
df = pd.read_csv("regress.results.csv")

if "${params.incl_unaligned}" != "false":
    df = df.query("feature != 'unaligned_reads'")

for param, param_df in df.groupby("parameter"):
    if param != "(Intercept)":
        (
            param_df
            .set_index("feature")
            .drop(columns=["parameter"])
            .to_csv(f"{param}.results.csv")
        )
    """
}

process plot_regress {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/association/", mode: 'copy', overwrite: true

    input:
    path "regress.results.csv"

    output:
    path "*"

    """#!/bin/bash
set -e
    
plot_regress.py
"""
}

process plot_metagenomes {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/association/${param}/", mode: 'copy', overwrite: true

    input:
    tuple val(param), path(stats), path(h5ad)

    output:
    path "*"
    path "${stats}"

    """#!/bin/bash
set -e
    
plot_metagenomes.py \
    --param "${param}" \
    --stats "${stats}" \
    --h5ad "${h5ad}" \
    --output_folder ./
"""
}