#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
Import require_param from 'lib/helpers.groovy'
Import help_message from 'lib/helpers.groovy'

// Import sub-workflow
include { deduplicate } from './modules/deduplicate'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    help_message(
        """
        Dedicated gene deduplication utility
        
        Clusters a collection of gene sequences by amino acid similarity and outputs
        the centroids of each cluster.

        Parameters:

        --genes                Folder(s) containing all gene FASTA files to deduplicate
        --output               Folder where output files will be written

        --cluster_similarity   Amino acid similarity used for clustering (ranges from 0.0 to 1.0)
                            (default: ${params.cluster_similarity})
        --cluster_coverage     Alignment coverage used for clustering (ranges from 0.0 to 1.0)
                            (default: ${params.cluster_coverage})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    require_param(params.output, "output")
    require_param(params.genes, "genes")

    // Get all of the files in the specified folder
    Channel
        .fromPath(params.genes)
        .set { gene_ch }

    // Run the deduplication utility on those files
    deduplicate(gene_ch)

}