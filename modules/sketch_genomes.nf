#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Find ORFs in a genome
process find_orfs {
    container "${params.container__emboss}"
    label "io_limited"

    input:
    tuple path(genome), val(uri)

    output:
    tuple path("${genome}.orfs.fasta"), val(uri)

    script:
    template "find_orfs.sh"

}

// Sketch the genome
process sketch {
    container "${params.container__mash}"
    label "mem_medium"
    publishDir "${uri}.msh", saveAs: {"${uri}.msh"}, mode: 'copy', overwrite: true

    input:
    tuple path(orfs), val(uri)

    output:
    path "*.msh"

    script:
    template "sketch.sh"

}
