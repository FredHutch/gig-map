#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { find_orfs; sketch } from './modules/sketch_genomes'
include { group_sketches as group_batches } from './modules/sketch_genomes' addParams(save_sketches: false)
include { group_sketches as group_all } from './modules/sketch_genomes' addParams(save_sketches: true)

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Sketch Genomes
        
        Create an easily searchable sketch of the amino acid sequences contained in the open reading
        frames of all genome files contained in a folder.

        Parameters:

        --genomes           Folder containing the set of genome nucleotide FASTAs to sketch
        --sketch_folder     Folder which should contain the combined sketches of all genomes (combined_genomes.msh)
        --recursive         Include files contained in subdirectories
                            (default: ${params.recursive})
        --minsize           Minimum size of open reading frames
                            (default: ${params.minsize})
        --gencode           Genetic code used for conceptual translation of genome sequences
                            (default: ${params.gencode})
        --sketch_size       Number of minmers to include in each sketch
                            (default: ${params.sketch_size})
        --kmer_size         Length of k-mers used for sketching
                            (default: ${params.kmer_size})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.genomes, "genomes")
    helpers.require_param(params.sketch_folder, "sketch_folder")

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
            "${path_base}.fa"
        ] )
        .ifEmpty { error "Cannot find any files in ${genome_folder}/" }
        .set { genomes_ch }

    // Find the open reading frames in a genome
    find_orfs(genomes_ch)

    // Make a sketch from each set of kmers
    sketch(find_orfs.out)

    // Group together the sketches in two rounds
    group_batches(
        sketch
            .out
            .collate(params.batchsize)
    )

    // Collect as a single queryable sketch
    group_all(
        group_batches
            .out
            .collect()
    )

    // Write out all of the file paths of the genomes
    // which were indexed
    genomes_ch
        .map { it -> "${it}" }
        .collectFile(
            newLine: true,
            name: "combined_genomes.txt",
            storeDir: "${params.sketch_folder}/"
        )

}