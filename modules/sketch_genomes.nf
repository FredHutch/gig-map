#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Find ORFs in a genome
process find_orfs {
    container "${params.container__emboss}"
    label "io_limited"

    input:
    path genome

    output:
    path "${genome}.orfs.fasta"

    script:
    template "find_orfs.sh"

}

// Sketch the genome
process sketch {
    container "${params.container__mash}"
    label "mem_medium"

    input:
    path orfs

    output:
    path "*.msh"

    script:
    template "sketch.sh"

}

// Calculate the similarity of a query to a genome
process search {
    container "${params.container__mash}"
    label "mem_medium"

    input:
    path query_msh
    path ref_msh

    output:
    path "*.tsv"

    script:
    template "search.sh"

}

// Reformat the results
process reformat {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.search_results}/", mode: 'copy', overwrite: true, pattern: "*.csv"

    input:
    path tsv

    output:
    path "*.csv"

    script:
    template "reformat.py"

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

// Combine a set of sketches into a single file
process group_sketches {
    container "${params.container__mash}"
    label "io_limited"
    publishDir "${params.sketch_folder}/", mode: 'copy', overwrite: true, enabled: params.save_sketches

    input:
    // Rename the inputs ordinally to prevent filename collisions
    path "inputs/?.msh"

    output:
    path "combined_genomes.msh"

    script:
    template "group_sketches.sh"

}