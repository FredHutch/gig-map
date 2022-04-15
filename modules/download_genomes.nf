#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the processes to run in this workflow
include {
    fetchFTP;
    parse_genome_csv;
    concatenate_annotations;
} from './processes/download' addParams(
    ftp_output_folder: "${params.output}/genomes",
    publishFTP: 'true',
)


workflow download_genomes {

    take:
    genome_manifests

    main:

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
            .toSortedList(),
        "genomes"
    )

    emit:
    genomes = fetchFTP.out
    annot = concatenate_annotations.out

}