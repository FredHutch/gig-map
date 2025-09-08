#!/usr/bin/env nextflow

// Entrypoint used to build trees from nucleotide gene sequences

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

include { build_trees } from './modules/build_trees'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Build phylogenetic trees from a collection of gene nucleotide sequences

        Parameters:

        --genes             Folder containing the set of gene amino acid FASTAs to deduplicate.
                            Every file must end with ".fasta.gz"
        --gene_bins         Grouping of genes into bins (e.g. gene_bins.csv)
        --output            Folder where output files will be written

        Notes:
            Each of the files in the --genes folder should be a FASTA file
            containing the nucleotide sequences of a single gene across multiple genomes.
            The file name (without the suffix) will be used as the gene identifier
            and should match the gene IDs in the --gene_bins file (column 'gene_id').

            The 'bin' column in the --gene_bins file will be used to group genes
            together for tree construction. Each unique bin will result in a single
            tree containing all genes assigned to that bin.
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.genes, "genes")
    helpers.require_param(params.gene_bins, "gene_bins")
    helpers.require_param(params.output, "output")

    // Make a channel with the input gene FASTA files
    Channel
        .fromPath("${params.genes}/*.fasta.gz")
        .ifEmpty { error "No gene FASTA files found in: ${params.genes}" }
        .set { gene_fastas }

    // Run the subworkflow to align the genes and build trees for each bin
    build_trees(
        gene_fastas,
        file(params.gene_bins, checkIfExists: true)
    )

}