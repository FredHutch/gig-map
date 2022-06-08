#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { sketch; search; reformat } from './modules/sketch_genomes' addParams(save_sketches: false)

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Search Genome Sketches
        
        Search one or more protein queries (multi-FASTAs) against all of the sketched
        genomes contained in a folder. It is expected that the genomes to search against
        will have been sketched using the 'sketch_genomes' utility in gig-map, which places
        an amino-acid *.msh file in the same folder as the genome it was generated from.

        Parameters:

        --query             Protein sequence(s) to search against those genomes, in multi-FASTA format.
                            Multiple queries may be specified with wildcard characters.
        --genome_sketches   File containing genome sketches to search against (*.msh)
        --search_results    Folder used to write results
        --sketch_size       Number of minmers to include in each query sketch
                            (default: ${params.sketch_size})
        --kmer_size         Length of k-mers used for sketching
                            (note: this must match the size used to sketch the genomes)
                            (default: ${params.kmer_size})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.query, "query")
    helpers.require_param(params.genome_sketches, "genome_sketches")
    helpers.require_param(params.search_results, "search_results")

    // Log the parameters being used
    log.info"""
    query:            ${params.query}
    genome_sketches:  ${params.genome_sketches}
    search_results:   ${params.search_results}
    sketch_size:      ${params.sketch_size}
    kmer_size:        ${params.kmer_size}
    """

    // Make a channel with the query or queries
    Channel
        .fromPath("${params.query}")
        .ifEmpty( { error "No inputs found at query location: ${params.query}" } )
        .set { query_ch }

    // Make a sketch from each set of the queries
    sketch(query_ch)

    // Count the number of minmers which overlap with all of the genomes in the query
    search(
        sketch.out,
        file("${params.genome_sketches}", checkIfExists: true)
    )

    // Reformat the results
    reformat(
        search.out
    )

}
