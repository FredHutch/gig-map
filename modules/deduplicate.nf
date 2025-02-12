#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import processes
include {
    cdhit_scatter;
    cdhit_gather;
    merge_cluster_membership;
    deduplicate_fasta_names;
    annotate_centroids;
    filter_genes;
    get_gene_annot;
} from './processes/deduplicate'

workflow deduplicate {

    take:
    genes_ch

    main:

    // Filter the genes by minimum amino acid length
    filter_genes(
        genes_ch
    )

    // Get the annotations for all of the genes
    get_gene_annot(
        genes_ch
    )

    // Deduplicate the gene sequences
    deduplicate_fasta_names(
        filter_genes
            .out
            .toSortedList()
    )

    // Run CD-HIT on all of the shards
    cdhit_scatter(
        deduplicate_fasta_names.out.flatten()
    )

    // Gather the CD-HIT results
    cdhit_gather(
        cdhit_scatter.out.fasta.toSortedList()
    )

    // Merge the cluster membership labels
    merge_cluster_membership(
        cdhit_scatter.out.membership.toSortedList(),
        cdhit_gather.out.membership
    )

    // Generate a simple annotation file for each centroid
    annotate_centroids(
        cdhit_gather.out.fasta,
        get_gene_annot
            .out
            .annot
            .toSortedList()
    )

    emit:
    fasta = cdhit_gather.out.fasta
    annot = annotate_centroids.out.annot

}