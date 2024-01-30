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
    parse_genome_csv_suffix: "_protein.faa.gz",
    ftp_output_folder: "${params.output}/genes",
    skip_missing_ftp: "true",
    publishFTP: 'true',
    ncbi_datasets_type: 'protein'
)

workflow download_genes {

    take:
    gene_manifests_csv
    gene_manifests_tsv
    
    main:

    // Read the contents of each manifest file
    parse_genome_csv(
        gene_manifests_csv
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
        gene_manifests_tsv
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
        "genes"
    )

    emit:
    genes = fetchFTP.out.mix(unpackDatasets.out.fasta)
    annot = concatenate_annotations.out
}
