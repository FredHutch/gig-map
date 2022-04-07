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

    output:
    path "*.html"

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
    ${params.render_options}
    """
}

workflow {

    // *****************
    // GENOME ALIGNMENTS
    // *****************

    // If the user specified custom genome alignments path
    if ( params.genome_alignments ){
        // Use that path for the annotations
        genome_alignments = params.genome_alignments
    
    // If the user did not specify a custom genome alignments path
    } else {
        // Use the default path in the project folder
        genome_alignments = "${params.project_folder}/genome_alignments/genomes.aln.csv.gz"
    }

    // Get that file
    Channel
        .fromPath(genome_alignments)
        .ifEmpty { error "No file found at ${genome_alignments}" }
        .set { genome_alignments }

    // ******************
    // GENOME ANNOTATIONS
    // ******************

    // If the user specified custom genome annotations path
    if ( params.genome_annotations ){
        // Use that path for the annotations
        genome_annotations = params.genome_annotations
    
    // If the user did not specify a custom genome annotations path
    } else {
        // Use the default path in the project folder
        genome_annotations = "${params.project_folder}/downloaded_genomes/genomes.annot.csv.gz"
    }

    // Get that file
    Channel
        .fromPath(genome_alignments)
        .ifEmpty { error "No file found at ${genome_alignments}" }
        .set { genome_alignments }

    // ****************
    // GENOME DISTANCES
    // ****************

    // If the user specified custom genome distmat path
    if ( params.genome_distmat ){
        // Use that path for the annotations
        genome_distmat = params.genome_distmat
    
    // If the user did not specify a custom genome distmat path
    } else {
        // Use the default path in the project folder
        genome_distmat = "${params.project_folder}/ani/distances.csv.gz"
    }

    // Get that file
    Channel
        .fromPath(genome_distmat)
        .ifEmpty { error "No file found at ${genome_distmat}" }
        .set { genome_distmat }

    // **********
    // GENE ORDER
    // **********

    // If the user specified custom gene order path
    if ( params.gene_order ){
        // Use that path for the annotations
        gene_order = params.gene_order
    
    // If the user did not specify a custom gene order path
    } else {
        // Use the default path in the project folder
        gene_order = "${params.project_folder}/genome_alignments/genomes.gene.order.txt.gz"
    }

    // Get that file
    Channel
        .fromPath(gene_order)
        .ifEmpty { error "No file found at ${gene_order}" }
        .set { gene_order }

    // ****************
    // GENE ANNOTATIONS
    // ****************

    // If the user specified custom gene annotations path
    if ( params.gene_annotations ){
        // Use that path for the annotations
        gene_annotations = params.gene_annotations
    
    // If the user did not specify a custom gene annotations path
    } else {
        // Use the default path in the project folder
        gene_annotations = "${params.project_folder}/deduplicated_genes/centroids.annot.csv.gz"
    }

    // Get that file
    Channel
        .fromPath(gene_annotations)
        .ifEmpty { error "No file found at ${gene_annotations}" }
        .set { gene_annotations }

    render(
        genome_alignments
        genome_annotations
        genome_distmat
        gene_order
        gene_annotations
    )

}