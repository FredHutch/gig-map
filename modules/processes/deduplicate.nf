// Deduplicate the names before running CDHIT
process deduplicate_fasta_names {
    container "${params.container__pandas}"
    label 'io_limited'
   
    input:
    path "input.genes.*.fasta.gz"
    
    output:
    path "deduplicated.genes.*.fasta.gz"
    
    script:
    """#!/bin/bash
deduplicate_fasta_names.py 'input.genes.*.fasta.gz' deduplicated.genes '${params.cluster_shard_size}'

    """
}

// Cluster genes by amino acid similarity, once per shard
process cdhit_scatter {
    container "${params.container__cdhit}"
    label 'mem_medium'
   
    input:
    path "deduplicated.genes.0.fasta.gz"
    
    output:
    path "centroids.faa.gz", emit: fasta
    path "centroids.membership.csv.gz", emit: membership
    
    script:
    template "cdhit.sh"
}

// Cluster genes by amino acid similarity
process cdhit_gather {
    container "${params.container__cdhit}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "*.faa.gz"
   
    input:
    path "deduplicated.genes.*.fasta.gz"
    
    output:
    path "centroids.faa.gz", emit: fasta
    path "centroids.membership.csv.gz", emit: membership
    
    script:
    template "cdhit.sh"
}

// Merge the membership information from all of the shards
process merge_cluster_membership {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path "scatter.membership.*.tsv.gz"
    path "gather.membership.tsv.gz"

    output:
    path "centroids.membership.csv.gz"

    script:
    """#!/bin/bash
merge_cluster_membership.py
    """
}


// Generate a simple annotation path for each centroid
process annotate_centroids {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path "clustered.genes.fasta.gz"
    path "input_annotations/"
    
    output:
    path "centroids.annot.csv.gz", emit: annot
    
    script:
    template "annotate_centroids.py"
}

// Get any annotations which are present in the FASTA headers
process get_gene_annot {
    container "${params.container__pandas}"
    label 'io_limited'
   
    input:
    path fasta
    
    output:
    path "${fasta}.annot.csv.gz", emit: annot
    
    script:
    template "get_gene_annot.py"
}

// Filter a set of genes by minimum amino acid length
process filter_genes {
    container "${params.container__pandas}"
    label 'io_limited'
    
    input:
        path input_fasta
    
    output:
        path "${input_fasta.name}.filtered.fasta.gz"
    
    script:
    template "filter_genes.py"

}