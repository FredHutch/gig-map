#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { find_orfs; sketch; collect } from './modules/sketch_genomes'

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
        .map { it -> [it, "${it}"]}
        .set { genomes_ch }

    // Find the open reading frames in a genome
    find_orfs(genomes_ch)

    // Make a sketch from each set of kmers
    sketch(find_orfs.out)

}