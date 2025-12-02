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

        --genes             Folder containing the set of gene nucleotide FASTAs to build trees from. 
                            This should be the align/genes/ folder output by build_pangenome.nf
        --gene_bins         Grouping of genes into bins (e.g. gene_bins.csv)
        --raxml_min_prop_sites
                            Filter out any positions which have gaps for this proportion of genomes
        --raxml_min_prop_genomes
                            Filter out any genomes which have gaps for this proportion of sites
        --min_raxml_genomes Only build alignments for genes with at least this number of filtered genomes
        --max_raxml_sites   Truncate the filtered MSA to contain no more than this number of sites

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