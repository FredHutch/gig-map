#!/usr/bin/env nextflow

// Default parameters are in `nextflow.config`

// Using DSL-2
nextflow.enable.dsl=2

// Import the various submodules used by gig-map
include { align_markers } from './align_markers'
include { ani } from './ani'
include { aggregate } from './aggregate'
include { render } from './render'

include { clean_genomes } from './processes/align_genomes'

workflow collect {

    take:

        genomes_ch
        genome_aln
        gene_order
        marker_genes


    main:

    // Clean up the genome formatting
    clean_genomes(
        genomes_ch
            .map {
                it -> [it.name, it]
            }
    )

    // Align the set of marker genes against these genomes,
    // combining the genes found commonly across genomes with
    // any genes that the user may have placed in markers/
    align_markers(
        clean_genomes.out,
        marker_genes
    )

    // Compute the ANI similarity of all genomes
    ani(
        clean_genomes.out
    )

    // Aggregate results
    aggregate(
        genome_aln,
        ani.out.distances,
        align_markers.out.distmat,
        gene_order
    )

    // Render the results as an interactive figure
    render(
        genome_aln,
        aggregate.out.genome_annot,
        ani.out.distances,
        gene_order,
        aggregate.out.gene_annot
    )

}
