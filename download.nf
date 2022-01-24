#!/usr/bin/env nextflow

// Default parameters are in `nextflow.config`

// Using DSL-2
nextflow.enable.dsl=2

// Import the various submodules used by gig-map
include { download_genomes } from './modules/download_genomes'
include { download_genes } from './modules/download_genes'


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/gig-map/download <ARGUMENTS>

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

}
