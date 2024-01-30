#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the processes to run in this workflow
include {
    fetchFTP;
    parse_genome_csv;
    parse_genome_tsv;
    getDatasets;
    unpackDatasets;
    concatenate_annotations;
} from './processes/download' addParams(
    ftp_output_folder: "${params.output}/genomes",
    publishFTP: 'true',
)


workflow download_genomes {

    take:
    genome_manifest_csv
    genome_manifest_tsv

    main:

    // Read the contents of each CSV manifest file
    parse_genome_csv(
        genome_manifest_csv
    )

    // Download each of the files
    fetchFTP(
        parse_genome_csv
            .out[0]
            .splitText()
            .map({it.trim()})
    )

    // Read the contents of each TSV manifest file
    parse_genome_tsv(
        genome_manifest_tsv
    )

    // Download each of the Datasets
    getDatasets(
        parse_genome_tsv
            .out
            .splitText()
            .map({it.trim()})
    )
    // Unpack the ZIP archive that is downloaded
    unpackDatasets(getDatasets.out)

    // Join the genome annotations from the NCBI table
    concatenate_annotations(
        parse_genome_csv
            .out[1]
            .toSortedList(),
        unpackDatasets
            .out
            .annot
            .toSortedList(),
        "genomes"
    )

    emit:
    genomes = fetchFTP.out.mix(unpackDatasets.out.fasta)
    annot = concatenate_annotations.out

}