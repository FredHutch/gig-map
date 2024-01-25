process bin_metagenomes {
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
    path "*"


    """#!/bin/bash
set -e

bin_metagenomes.py \
    --read_alignments "${read_alignments}" \
    --gene_bins "${gene_bins}" \
    --genome_groups "${genome_groups}" \
    --group_profile "${group_profile}" \
    --metadata "${metadata}" \
    --category ${params.category} \
    --min_n_reads ${params.min_n_reads} \
    --min_n_genes ${params.min_n_genes} \
    --output_folder ./
    """
}