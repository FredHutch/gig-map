#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Make sure that the manifest is valid
process parse_manifest {
    container "${params.container__pandas}"
    label "io_limited"

    input:
    path "raw.manifest.csv"
    path "read_alignments.csv.gz"

    output:
    path "parsed.manifest.csv"

    script:
    template "parse_manifest.py"

}

// Break up the gene abundances
process shard_genes {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    file "manifest.csv"
    file "read_alignments.csv.gz"

    output:
    file "readcounts.*.csv.gz"

    script:
    template "shard_genes.py"

}

// Test for differences between samples
process regress {
    container "${params.container__regress}"
    label "cpu_high"
    
    input:
    file metadata_csv
    file readcounts_csv_gz

    output:
    path "regress.results.csv", optional: true

    script:
    template "regress.Rscript"

}

// Join the results
process join {
    container "${params.container__pandas}"
    publishDir "${params.output}", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    file "regress.results.*.csv"

    output:
    file "*.csv.gz"

    script:
    template "join.py"

}