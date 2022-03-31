#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the processes to run in this workflow
include {
    fetchFTP;
    parse_genome_csv;
    concatenate_annotations;
} from './processes/download' addParams(
    parse_genome_csv_suffix: "_protein.faa.gz",
    ftp_output_folder: params.output,
    skip_missing_ftp: "true",
    publishFTP: 'true',
)

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

// Set the default parameters
params.genome_csv = false
params.output = false

workflow download_genes {

    take:
    gene_manifests
    
    main:

    // Read the contents of each manifest file
    parse_genome_csv(
        gene_manifests
    )

    // Download each of the files
    fetchFTP(
        parse_genome_csv
            .out[0]
            .splitText()
            .map({it.trim()})
    )

    // Join the genome annotations from the NCBI table
    concatenate_annotations(
        parse_genome_csv
            .out[1]
            .toSortedList(),
        "genes"
    )

    emit:
    genes = fetchFTP.out
}

workflow {

    // If the --genome_csv or --output flags were not set
    if ( params.output == false or params.genome_csv == false ){

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