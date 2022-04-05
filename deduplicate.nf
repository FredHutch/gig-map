#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflow
include { deduplicate } from './modules/deduplicate'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Dedicated gene deduplication utility
        
        Clusters a collection of gene sequences by amino acid similarity and outputs
        the centroids of each cluster.

        Parameters:

        --genes                Wildcard path indicating the amino acid gene FASTA files to analyze
        --output               Folder where output files will be written

        --cluster_similarity   Amino acid similarity used for clustering (ranges from 0.0 to 1.0)
                               (default: ${params.cluster_similarity})
        --cluster_coverage     Alignment coverage used for clustering (ranges from 0.0 to 1.0)
                               (default: ${params.cluster_coverage})
        --min_gene_length      Minimum amino acid length threshold used to filter genes
                               (default: ${params.min_gene_length})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genes, "genes")

    // Get all of the files in the specified folder
    Channel
        .fromPath(params.genes)
        .ifEmpty { error "Cannot find any files at '${params.genes}'" }
        .set { gene_ch }

    // Run the deduplication utility on those files
    deduplicate(gene_ch)

}