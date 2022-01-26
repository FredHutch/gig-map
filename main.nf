#!/usr/bin/env nextflow

// Default parameters are in `nextflow.config`

// Using DSL-2
nextflow.enable.dsl=2

// Import the various submodules used by gig-map
include { download_genomes } from './modules/download_genomes'
include { download_genes } from './modules/download_genes'
include { deduplicate } from './modules/deduplicate'
include { align_genomes } from './modules/align_genomes'
include { align_markers } from './modules/align_markers'
include { ani } from './modules/ani'
include { aggregate } from './modules/aggregate'
include { render } from './modules/render'


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/gig-map <ARGUMENTS>

    For detailed instructions, please visit https://github.com/FredHutch/gig-map/
    or the Readme.md file provided in the downloaded repository.

    """.stripIndent()
}


workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // Download genomes and genes specified by the user
    // in the genome_tables/ and gene_tables/ folders
    // (or the corresponding params)
    download_genomes()
    download_genes()

    // Deduplicate all input genes to get centroids for alignment
    deduplicate(
        download_genes
            .out
            .genes
    )

    // Align the genomes with the genes
    align_genomes(
        download_genomes
            .out
            .genomes,
        deduplicate
            .out
            .fasta
    )

    // Align the set of marker genes against these genomes,
    // combining the genes found commonly across genomes with
    // any genes that the user may have placed in markers/
    align_markers(
        align_genomes
            .out
            .clean_genomes,
        align_genomes
            .out
            .markers
    )

    // Compute the ANI similarity of all genomes
    ani(
        align_genomes
            .out
            .clean_genomes
    )

    // Aggregate results
    aggregate(
        align_genomes.out.concat_alignments,
        ani.out.distances,
        align_markers.out.distmat,
        align_genomes.out.gene_order
    )

    // // Render the results as an interactive figure
    // render(
    //     align_genomes.out.concat_alignments,
    //     download_genomes.out.annot,
    //     ani.out.distances,
    //     align_genomes.out.gene_order,
    //     deduplicate.out.annot
    // )

}
