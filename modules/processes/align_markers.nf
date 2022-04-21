// Make the multiple sequence alignment
process combine_markers {
    container "${params.container__clustal}"
    label "mem_medium"
    publishDir "${params.output}/msa/", mode: 'copy', overwrite: true, pattern: "*.msa.gz"
    publishDir "${params.output}/distmat/", mode: 'copy', overwrite: true, pattern: "*.distmat"

    input:
        path unaligned_fasta

    output:
        path "*.distmat", emit: distmat
        path "*.msa.gz", emit: msa

"""#!/bin/bash

set -e

# Make a local copy of the decompressed FASTA
gunzip -c ${unaligned_fasta} > input.fasta

# Run Clustal-Omega
clustalo \
    --in input.fasta \
    -t DNA \
    --full \
    --distmat-out="${unaligned_fasta.name.replaceAll('.fasta.gz', '')}.distmat" \
    --out="${unaligned_fasta.name.replaceAll('.fasta.gz', '')}.msa" \
    --threads ${task.cpus} \
    --verbose \
    --outfmt=fasta

gzip *msa

"""
}

// Build an ML tree from the MSA
process raxml {
    container "${params.container__raxml}"
    label 'mem_medium'
    publishDir "${params.output}/raxml/", mode: 'copy', overwrite: true

    input:
        path aln_fasta
    
    output:
        path "${aln_fasta.name.replaceAll('.gz', '')}.raxml.bestTree"
        path "${aln_fasta.name.replaceAll('.gz', '')}.raxml.bestModel"
        path "${aln_fasta.name.replaceAll('.gz', '')}.raxml.log"
    
    script:
    template 'raxml.sh'
}