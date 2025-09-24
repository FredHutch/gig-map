#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { bin_metagenomes } from './modules/bin_metagenomes'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Summarize a set of metagenome alignments using a set of binned genes.
        Gene binning can be performed using co-occurrence across genomes, and
        so the binned genes represent genetic elements which are observed in a
        correlated manner across a pangenome collection.
        By applying the results of gene binning to the metagenome alignment
        results using the same gene catalog, the diversity of a collection of
        organisms can be collapsed into units of those co-occurrent genetic
        elements.

        Parameters:

        --read_alignments   Alignments from align_reads.nf (e.g. read_alignments.csv.gz)
        --gene_bins         Grouping of genes into bins from bin_genes.nf (e.g. gene_bins.csv)
        --group_profile     Gene content of genome groups from bin_genes.nf (e.g. group_profile.csv)
        --genome_groups     Grouping of genomes into groups (e.g. genome_groups.csv)
        --centroids_length  Length of each centroids sequence (e.g. centroids_length.csv.gz)
        --output            Folder where output files will be written

        --metadata          Optional: Metadata table used to compare samples (CSV)
        --formula           Optional: Column(s) from metadata table used for comparison
        --incl_unaligned    Include unaligned reads in the comparison (default: false)

        --min_n_reads       Exclude any samples with fewer than this number of reads aligned (default: 0)
        --min_n_genes       Exclude any samples with fewer than this number of genes detected (default: 0)

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.read_alignments, "read_alignments")
    helpers.require_param(params.gene_bins, "gene_bins")
    helpers.require_param(params.group_profile, "group_profile")
    helpers.require_param(params.centroids_length, "centroids_length")
    helpers.require_param(params.output, "output")

    // Get the read alignments
    read_alignments = file(params.read_alignments, checkIfExists: true)

    // Get the gene bins
    gene_bins = file(params.gene_bins, checkIfExists: true)

    // Get the genome groups
    genome_groups = file(params.genome_groups, checkIfExists: true)

    // Get the centroid length CSV
    centroids_length = file(params.centroids_length, checkIfExists: true)

    // Get the genome groups
    group_profile = file(params.group_profile, checkIfExists: true)

    // Get the metadata table
    metadata = file(params.metadata, checkIfExists: true)

    // Bin the metagenomes
    bin_metagenomes(
        read_alignments,
        centroids_length,
        gene_bins,
        genome_groups,
        group_profile,
        metadata
    )

}