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
    publishDir "${uri}.msh", saveAs: {"${uri}.msh"}, mode: 'copy', overwrite: true, enabled: params.save_sketches

    input:
    tuple path(orfs), val(uri)

    output:
    tuple path("*.msh"), val(uri)

    script:
    template "sketch.sh"

}

// Calculate the similarity of a query to a genome
process search {
    container "${params.container__mash}"
    label "mem_medium"

    input:
    tuple path(query_msh), val(query_filename), val(genome_uri), path(genome), path(genome_msh)

    output:
    tuple val(query_filename), path("*.tsv")

    script:
    template "search.sh"

}

// Collect the results
process collect {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.search_results}/", mode: 'copy', overwrite: true, enabled: true, pattern: "*.csv"

    input:
    // Rename the inputs ordinally to prevent filename collisions
    tuple val(query_filename), path("inputs/*.tsv")

    output:
    path "${query_filename}.csv", emit: csv
    tuple val(query_filename), path("overlapping.txt"), optional: true, emit: overlapping

    script:
    template "collect.py"

}

// Save the genome files which match a particular query
process save_genomes {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.search_results}/${query_filename}/", mode: 'copy', overwrite: true

    input:
    // Rename the inputs ordinally to prevent filename collisions
    tuple val(query_filename), path("matching_genomes/")

    output:
    path "matching_genomes/*", includeInputs: true

    """#!/bin/bash

    set -e
    ls -lahtr matching_genomes/
    """

}