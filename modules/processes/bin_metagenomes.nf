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

    script:
    template "bin_summary.py"

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

// Split up the corncob results as inputs for plotting
process split {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
    path "corncob.results.csv"

    output:
    path "*.results.csv"

    """#!/usr/bin/env python3
import pandas as pd
df = pd.read_csv("corncob.results.csv")

if "${params.incl_unaligned}" != "false":
    df = df.query("feature != 'unaligned_reads'")

for param, param_df in df.groupby("parameter"):
    if not param.startswith("mu.") or param.endswith("(Intercept)"):
        continue
    
    param = param[3:]
    (
        param_df
        .set_index("feature")
        .drop(columns=["parameter"])
        .to_csv(f"{param}.results.csv")
    )
    """
}

process plot {
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