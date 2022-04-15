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
    ftp_output_folder: "${params.output}/genes",
    skip_missing_ftp: "true",
    publishFTP: 'true',
)

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
