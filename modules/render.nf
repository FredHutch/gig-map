#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Define the process for rendering an HTML file from an RDB input
process render {
    container "${params.container__gigmap}"
    memory "${params.render_mem_gbs}.GB"
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path genome_alignments
    path genome_annotations
    path genome_distmat
    path gene_order
    path gene_annotations
    path render_options

    output:
    path "gigmap*"

    script:
    """#!/bin/bash

set -e

# Get the timestamp used for the file
NOW=\$(date +"%Y-%m-%d-%H-%M-%S")

gig-map-cli \
    --genomeTree-distmat "${genome_distmat}" \
    --genomeHeatmap-csv "${genome_alignments}" \
    --genomeAnnot-csv "${genome_annotations}" \
    --geneAnnot-csv "${gene_annotations}" \
    --geneAnnot-order "${gene_order}" \
    --output-folder ./ \
    --output-prefix gigmap.\$NOW \
    \$(cat render_options | jq '')

cp .command.sh gigmap.\$NOW.sh
cp .command.log gigmap.\$NOW.log
    """
}