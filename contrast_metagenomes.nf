#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { contrast_metagenomes } from './modules/contrast_metagenomes'
include { centroids_length } from './modules/processes/align_reads'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Perform statistical analysis comparing the proportional abundance of bins between
        contrasting groups of metagenomes.
        This entrypoint should be used on the outputs of the "align_pangenome" workflow.

        The --formula parameter should reference annotations available in the --metadata
        CSV, where the first column includes sample IDs that correspond to those used
        when aligning the FASTQ inputs against the pangenome.

        Parameters:

        --read_alignments   Alignments from align_reads.nf (e.g. read_alignments.csv.gz)
        --gene_bins         Grouping of genes into bins from bin_genes.nf (e.g. gene_bins.csv)
        --centroids_length  Length of each centroids sequence (e.g. centroids_length.csv.gz)

        --metadata          Metadata table used to compare samples (CSV)
        --formula           Column(s) from metadata table used for comparison (multiple columns specified like "colA + colB")

        --output            Folder where output files will be written

        --min_n_reads       Exclude any samples with fewer than this number of reads aligned (default: 0)
        --min_n_genes       Exclude any samples with fewer than this number of genes detected (default: 0)

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.read_alignments, "read_alignments")
    helpers.require_param(params.gene_bins, "gene_bins")
    helpers.require_param(params.metadata, "metadata")
    helpers.require_param(params.formula, "formula")
    helpers.require_param(params.output, "output")

    // Get the read alignments file(s)
    Channel
        .fromPath( params.read_alignments.split(',').toList() )
        .toSortedList()
        .set { read_alignments }

    // Get the gene bins
    gene_bins = file(params.gene_bins, checkIfExists: true)

    // Get the centroids lengths
    centroids_length_file = file(params.centroids_length)

    // If the centroids lengths do not exist, build that file
    if(centroids_length_file.exists() == false){
        centroids_faa = file(params.centroids_faa, checkIfExists: true)
        centroids_length(centroids_faa)
        centroids_length_file = centroids_length.out
    }

    // Get the metadata table
    metadata = file(params.metadata, checkIfExists: true)

    // Bin the metagenomes
    contrast_metagenomes(
        read_alignments,
        centroids_length_file,
        gene_bins,
        metadata
    )

}