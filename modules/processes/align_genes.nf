// Filter the gene MSAs based on minimum thresholds for site and genome presence
process filter_msas {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
        path "input/"
    output:
        path "*.msa.gz", emit: filtered_msa, optional: true
    script:
        template 'filter_msas.py'
}

// Make the multiple sequence alignment
process gene_msa {
    container "${params.container__clustal}"
    label "mem_medium"
    publishDir "${params.output}/msa/", mode: 'copy', overwrite: true, pattern: "*.msa.gz"

    input:
        path unaligned_fasta

    output:
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
    --out="${unaligned_fasta.name.replaceAll('.fasta.gz', '')}.msa" \
    --threads ${task.cpus} \
    --verbose \
    --outfmt=fasta

gzip *msa

"""
}


process calc_distmat {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/distmat/", mode: 'copy', overwrite: true
    input:
        path msa_fasta
    output:
        path "${msa_fasta.name.replaceAll('.msa.gz', '')}.distmat.csv.gz", emit: distmat
    script:
    template 'calc_distmat.py'
}


process merge_bins {
    // Combine the distance matrices for each bin
    // into a single distance matrix per bin
    // Also, combine the MSAs for each bin into a single MSA per bin
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/bin_distmat/", mode: 'copy', overwrite: true, pattern: "*.distmat.csv.gz"
    input:
        path "msas/"
        path "distmats/"
        path gene_bins
    output:
        path "*.distmat.csv.gz", emit: distmat
        path "*.msa.gz", emit: msa
    script:
    template 'merge_bins.py'
}


// Build an ML tree from the MSA
process raxml {
    container "${params.container__raxml}"
    label 'cpu_high'
    publishDir "${params.output}/raxml/", mode: 'copy', overwrite: true

    input:
        path aln_fasta
    
    output:
        path "${aln_fasta.name.replaceAll('.gz', '')}.raxml.*"
    
    script:
    template 'raxml.sh'
}