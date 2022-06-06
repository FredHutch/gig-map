#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { sketch; search; collect; save_genomes } from './modules/sketch_genomes' addParams(save_sketches: false)

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
        --genomes           Folder containing the set of sketched genomes
        --search_results    Folder used to write results
        --save_overlapping  If specified, copy all genome files which have at least this number of overlapping minmers
        --recursive         Include genome sketches contained in subdirectories within the --genomes folder
                            (default: ${params.recursive})
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
    helpers.require_param(params.genomes, "genomes")
    helpers.require_param(params.search_results, "search_results")

    // Log the parameters being used
    log.info"""
    query:            ${params.query}
    genomes:          ${params.genomes}
    search_results:   ${params.search_results}
    save_overlapping: ${params.save_overlapping}
    recursive:        ${params.recursive}
    sketch_size:      ${params.sketch_size}
    kmer_size:        ${params.kmer_size}
    """

    // Remove any trailing slash from the genome folder
    def genome_folder = params.genomes.replaceAll('/$', '')

    // Use * or ** depending on the --recursive flag
    def path_base = ("${params.recursive}" == true) ? "${genome_folder}/**" : "${genome_folder}/*"

    // Get all of the genomes
    Channel
        // Check all of the possible file extensions
        .fromPath( [
            "${path_base}.fna.gz",
            "${path_base}.fasta.gz",
            "${path_base}.fa.gz",
            "${path_base}.fna",
            "${path_base}.fasta",
            "${path_base}.fa",
            "${path_base}.msh"
        ] )
        .ifEmpty { error "Cannot find any files in ${genome_folder}/" }
        // Key each file to the version without a trailing .msh
        .map { it -> ["${it}".replaceAll(/.msh$/, ''), it] }
        // Pair each file with its prefix
        .groupTuple()
        // Omit any files which do not have a sketch
        .filter { it.size() == 2 }
        // Flatten to a tuple of [ genome_uri, genome_file, sketch_file ]
        .map { it -> [it[0], it[1][0], it[1][1]] }
        .set { genomes_ch }

    // Make a channel with the query or queries
    Channel
        .fromPath("${params.query}")
        .ifEmpty( { error "No inputs found at query location: ${params.query}" } )
        // Keep track of the file alongside its name
        .map { it -> [it, "${it.name}"] }
        .set { query_ch }

    // Make a sketch from each set of the queries
    sketch(query_ch)

    // Search the queries against the genomes
    sketch
        .out
        .combine(genomes_ch)
        .set { search_ch }

    // The content of search_ch is:
        // query .msh
        // query filename
        // genome uri
        // genome file
        // genome .msh

    // Count the number of minmers which overlap
    search(search_ch)

    // Collect the results
    collect(
        search
            .out
            .groupTuple()
    )

    // If the user has asked to copy out all genomes with overlapping minmers
    if ( "${params.save_overlapping}" != "false" ){
        save_genomes(
            collect
                .out
                .overlapping
                .map {
                    it -> [
                        it[0],
                        it[1].splitText()
                    ]}
                .transpose()
                // Get the files which match for each query
                .map { it -> [it[0], file(it[1].trim(), checkIfExists: true)] }
                // Save those genomes to a folder named for the query
                .groupTuple()
        )
    }
    
}
