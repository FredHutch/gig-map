#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import the subworkflow to run
include { download_genomes } from './modules/download_genomes'

// Standalone entrypoint
workflow {

    helpers.help_message(
        """
        Download genomes from a collection of genomes in the NCBI Genome database

        Downloads the genome sequences in nucleotide FASTA format from a collection
        of genomes in the NCBI Genome database.
        
        The legacy format for a table of genomes to be downloaded from NCBI
        can be accessed in CSV format from this website:
        https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/
        Specify inputs from this website as --genome_csv

        Alternatively, the updated format for NCBI genomes using the Datasets
        API can be accessed in TSV format from this website:
        https://www.ncbi.nlm.nih.gov/datasets/genome/
        Specify inputs from this website as --genome_tsv

        Parameters:

        --genome_csv           Table of genomes to be downloaded (CSV)
        --genome_tsv           Table of genomes to be downloaded (TSV)
        --output               Folder where output files will be written

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")

    // Parse the CSV path(s), if provided
    Channel
        .fromPath( params.genome_csv.split(',').toList() )
        .set { genome_manifest_csv }

    // Parse the TSV path(s), if provided
    Channel
        .fromPath( params.genome_tsv.split(',').toList() )
        .set { genome_manifest_tsv }

    // Download the genomes for each genome
    download_genomes(
        genome_manifest_csv,
        genome_manifest_tsv
    )

}