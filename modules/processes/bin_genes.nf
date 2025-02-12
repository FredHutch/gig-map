process bin_genes {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path genome_aln
    path gene_annot
    path genome_annot

    output:
    path "*"

    """#!/bin/bash
set -e

bin_genes.py \
    --genome_aln "${genome_aln}" \
    --gene_annot "${gene_annot}" \
    --genome_annot "${genome_annot}" \
    --min_coverage "${params.min_coverage}" \
    --min_identity "${params.min_identity}" \
    --min_genomes_per_gene "${params.min_genomes_per_gene}" \
    --max_dist_genes "${params.max_dist_genes}" \
    --min_bin_size "${params.min_bin_size}" \
    --max_dist_genomes "${params.max_dist_genomes}" \
    --gene_proximity_enabled "${params.gene_proximity_enabled}" \
    --gene_proximity_threshold "${params.gene_proximity_threshold}" \
    """
}