#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import the subworkflow to run
include { download_genes } from './modules/download_genes'

// Function which prints help message text
def helpMessage() {
    log.info.stripIndent()
}


workflow {

    helpers.help_message(
        """
        Download genes from a collection of genomes in the NCBI Genome database
        
        Downloads the gene sequences in amino acid FASTA format from a collection
        of genomes in the NCBI Genome database. A table of genomes to be downloaded
        can be accessed from the website:
        https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/

        Parameters:

        --genome_csv           Table of genomes to be downloaded (CSV)
        --output               Folder where output files will be written

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genome_csv, "genome_csv")

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