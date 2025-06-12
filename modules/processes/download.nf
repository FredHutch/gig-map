// Fetch a file via FTP
process fetchFTP {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.ftp_output_folder}", mode: 'copy', overwrite: true, enabled: "${params.publishFTP}" == "true"

    maxForks params.ftp_threads

    input:
        val ftp_url
    
    output:
        path "*" optional params.skip_missing_ftp == "true"

    script:
    template "fetchFTP.py"
    
}


// Parse the NCBI Genome Browser CSV 
process parse_genome_csv {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        path "input.csv"
    
    output:
        path "url_list.txt"
        path "genome_annotations.csv.gz"

    script:
    template "parse_genome_csv.py"

}

// Parse the NCBI Genome Browser TSV 
process parse_genome_tsv {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        path "input.tsv"
    
    output:
        path "acc_list.txt"

    script:
    template "parse_genome_tsv.py"

}

process getDatasets {
    container "${params.container__datasets}"
    label "io_limited"
    tag "${dataset_acc}"

    maxForks params.ftp_threads

    input:
        val dataset_acc

    output:
        tuple val(dataset_acc), path("ncbi_dataset.zip"), optional: true

    script:
    template "getDatasets.sh"

}

process unpackDatasets {
    container "${params.container__famli}"
    label "io_limited"
    publishDir "${params.ftp_output_folder}", mode: 'copy', overwrite: true, enabled: "${params.publishFTP}" == "true", pattern: "*.f*"

    maxForks params.ftp_threads

    input:
        tuple val(dataset_acc), path("ncbi_dataset.zip")

    output:
        path "*.f*", emit: fasta
        path "assembly_data_report.jsonl", emit: annot

    """#!/bin/bash
set -e
unzip ncbi_dataset.zip
mv ncbi_dataset/data/assembly_data_report.jsonl ./
mv ncbi_dataset/data/*/* ./
[ -s protein.faa ] && mv protein.faa ${dataset_acc}_protein.faa
gzip *.f*
    """

}


// Combine all of the genome annotations into a single file
process concatenate_annotations {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path "annotations/*.csv.gz"
    path "annotations/*.jsonl"
    val output_prefix

    output:
    path "${output_prefix}.annot.csv.gz"

    script:
    template "concatenate_annotations.py"
}