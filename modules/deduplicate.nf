#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import processes
include {
    cdhit;
    annotate_centroids;
    filter_genes;
} from './processes/deduplicate'

workflow deduplicate {

    take:
    genes_ch

    main:

    // Filter the genes by minimum amino acid length
    filter_genes(
        genes_ch
    )


    // Run CD-HIT on all of the gene sequences
    cdhit(
        filter_genes
            .out
            .toSortedList()
    )

    // Generate a simple annotation file for each centroid
    annotate_centroids(
        cdhit.out.fasta
    )

    emit:
    fasta = cdhit.out.fasta
    annot = annotate_centroids.out.annot

}