#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the processes to run in this workflow
include {
    fetchFTP;
    parse_genome_csv;
    concatenate_annotations;
} from './processes' addParams(
    ftp_output_folder: "${params.project_folder}/downloaded_genomes",
    publishFTP: 'true',
)


workflow download_genomes {

    // If the --genome_tables flag is set, use that path
    if ( params.genome_tables ){
        genome_tables = params.genome_tables
    } else {
        // Otherwise use the default path in the project folder
        genome_tables = "${params.project_folder}/genome_tables/*.csv"
    }

    // Get the CSV files
    Channel
        .fromPath( genome_tables )
        .ifEmpty { error "Cannot find any files matching the wildcard '${genome_tables}'" }
        .set { genome_manifests }

    // Read the contents of each manifest file
    parse_genome_csv(
        genome_manifests
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
    genomes = fetchFTP.out

}