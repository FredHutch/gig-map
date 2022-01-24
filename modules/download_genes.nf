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
    ftp_output_folder: "${params.project_folder}/downloaded_genes",
    skip_missing_ftp: "true",
    publishFTP: 'true',
)


workflow download_genes {

    // If the --gene_tables flag is set, use that path
    if ( params.gene_tables ){
        gene_tables = params.gene_tables
    } else {
        // Otherwise use the default path in the project folder
        gene_tables = "${params.project_folder}/gene_tables/*.csv"
    }

    // Get the CSV files
    Channel
        .fromPath( gene_tables )
        .ifEmpty { error "Cannot find any files matching the wildcard '${gene_tables}'" }
        .set { gene_manifests }

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
            .toSortedList()
    )

    emit:
    genes = fetchFTP.out
}