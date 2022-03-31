#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the subworkflow to run
include { download_genes } from '.modules/download_genes'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Download genes from a collection of genomes in the NCBI Genome database
    
    Downloads the gene sequences in amino acid FASTA format from a collection
    of genomes in the NCBI Genome database. A table of genomes to be downloaded
    can be accessed from the website:
    https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/

    Parameters:

    --genome_csv           Table of genomes to be downloaded (CSV)
    --output               Folder where output files will be written

    """.stripIndent()
}

workflow {

    // If the --genome_csv or --output flags were not set
    if ( params.output == false || params.genome_csv == false ){

        // Print the help message
        helpMessage()

        // Raise an error
        exit 1

    } else {

        // Parse the CSV and raise an error if the path is not valid
        Channel
            .fromPath( params.genome_csv )
            .ifEmpty { error "Cannot find any file at '${params.genome_csv}'" }
            .set { gene_manifest }

        // Download the genes for each genome
        download_genes(
            gene_manifest
        )

    }

}