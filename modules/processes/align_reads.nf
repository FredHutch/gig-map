// Join pairs of FASTA files
process join_read_pairs {
    container "${params.container__gigmap}"
    label 'io_limited'
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path("inputs/")
    
    output:
    tuple val(sample_name), path("${sample_name}.${params.reads_suffix}")

    """#!/bin/bash
set -e
cat inputs/* > "${sample_name}.${params.reads_suffix}"
    """
}

// Count the number of reads per specimen
process count_reads {
    container "${params.container__pandas}"
    label 'io_limited'
    tag "${sample_name}"
    
    input:
    // Place all input files in an input/ folder, naming with a simple numeric index
    tuple val(sample_name), path("input/input*.fastq.gz")
    
    output:
    tuple val(sample_name), path("${sample_name}.num_reads.txt")

    """#!/bin/bash

# Decompress all FASTQ files
# Concatenate them
# Count the number of lines, divided by 4
gunzip -c input/input*.fastq.gz | awk 'NR % 4 == 1' | wc -l > ${sample_name}.num_reads.txt
    """
}

// Concatenate the read counts from multiple lanes
process concat_readcounts {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
    tuple val(sample_name), path("inputs/*.num_reads.txt")

    output:
    path "${sample_name}.num_reads.txt"

    script:
    template "concat_readcounts.py"
    
}

// Align short reads against a gene catalog in amino acid space
process diamond {
    container "${params.container__diamond}"
    label 'mem_medium'
    tag "${sample_name}"
    
    input:
    // Place all input files in an input/ folder, naming with a simple numeric index
    file refdb
    tuple val(sample_name), path("input/input*.fastq.gz")
    
    output:
    tuple val(sample_name), path("${sample_name}.aln.gz"), emit: aln
    path "*.log", emit: log

    shell:
    template "align_reads.sh"

}

process concat_aln {
    container "${params.container__pandas}"

    input:
    tuple val(sample_name), path("inputs/*.aln.gz")

    output:
    tuple val(sample_name), path("${sample_name}.aln.gz")

    script:
    """cat inputs/* > ${sample_name}.aln.gz"""
}

// Group together a collection of DIAMOND logs
process diamond_logs {
    container "${params.container__pandas}"
    label 'io_limited'
   
    input:
    path "*"

    output:
    path "alignment_logs.csv"

    script:
    template "diamond_logs.py"

}

// Filter out any samples which have an insufficient number of alignments
process filter_aln {
    container "${params.container__pandas}"
    label 'io_limited'
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(aln)
    
    output:
    tuple val(sample_name), path("${aln}"), optional: true

    script:
    template "filter_aln.py"

}

// Filter the alignments with the FAMLI algorithm
process famli {
    container "${params.container__famli}"
    label 'mem_medium'
    publishDir "${params.output}/alignments/", mode: "copy", overwrite: true
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), file(input_aln)
    
    output:
    path "${sample_name}.json.gz", emit: json
    path "${sample_name}.log", emit: log

    script:
    template "famli.sh"
}

// Summarize the amount of time used for FAMLI
process famli_logs {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/logs/", mode: "copy", overwrite: true
    
    input:
    path "*"
    
    output:
    path "famli_logs.csv"

    script:
    template "famli_logs.py"
}


// Combine the outputs from short read alignment across all specimens
process gather {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "famli/*"
    path "read_counts/*"
    
    output:
    path "read_alignments.csv.gz"

    """
    gather_alignments.py
    """   
}
