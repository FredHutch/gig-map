include {
    cdhit;
    annotate_centroids;
} from './processes/deduplicate'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Dedicated gene deduplication utility
    
    Clusters a collection of gene sequences by amino acid similarity and outputs
    the centroids of each cluster.

    """.stripIndent()
}


workflow deduplicate {

    take:
    // Genomes which have been previously downloaded
    downloaded_genes

    main:

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    

    // The user must specify an output folder
    if (!params.project_folder){
        log.info"""

        -----------------------
        MISSING --project_folder
        -----------------------

        """.stripIndent()
        helpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // If the --genes param has been set by the user
    if ( params.genes ) {

        // Get all of the files which are specified
        Channel
            .fromPath(params.genes)
            .set { gene_ch }

    // Otherwise, if the --genes param has not been set
    } else {

        // Read files from the genes/ subfolder in the project folder
        Channel
            .fromPath("${params.project_folder}/genes/*.gz")
            .set { gene_ch }

    }

    // Combine those files with the previously-downloaded genes (if any)
    gene_ch
        .mix(downloaded_genes)
        .ifEmpty { error "Cannot find any genes for deduplication" }
        .set { all_gene_ch }

    // Run CD-HIT on all of the gene sequences
    cdhit(
        all_gene_ch.toSortedList()
    )

    // Generate a simple annotation file for each centroid
    annotate_centroids(
        cdhit.out.fasta
    )

    emit:
    fasta = cdhit.out.fasta
    annot = annotate_centroids.out.annot

}