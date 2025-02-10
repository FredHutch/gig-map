// Make the multiple sequence alignment
process combine_markers {
    container "${params.container__clustal}"
    label "mem_medium"
    publishDir "${params.output}/msa/", mode: 'copy', overwrite: true, pattern: "*.msa.gz"
    publishDir "${params.output}/distmat/", mode: 'copy', overwrite: true, pattern: "*.distmat"

    input:
        path unaligned_fasta

    output:
        path "*.distmat", emit: distmat, optional: true
        path "*.msa.gz", emit: msa, optional: true

"""#!/bin/bash

set -e

# Make a local copy of the decompressed FASTA
gunzip -c ${unaligned_fasta} > input.fasta

# If there are fewer than 4 sequences, skip the alignment
if [ \$(grep -c ">" input.fasta) -lt 4 ]; then
    echo "Fewer than 4 sequences in the input, skipping alignment"
    exit 0
fi

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