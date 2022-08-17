#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { ani } from './modules/ani'

include { clean_genomes } from './modules/processes/align_genomes'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Compute ANI for a set of genomes
        
        Parameters:

        --genomes           Folder containing the set of genome nucleotide FASTAs to compare
        --output            Folder where output files will be written
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genomes, "genomes")

    // Remove any trailing slash from the genome folder
    genome_folder = params.genomes.replaceAll('/$', '')

    // Get all of the genomes
    Channel
        .fromPath("${genome_folder}/*")
        .ifEmpty { error "Cannot find any files at ${genome_folder}/*" }
        .set { genomes_ch }

    // Clean all of those genomes
    clean_genomes(
        genomes_ch
            .map {
                it -> [it.name, it]
            }
    )

    // Calculate the ANI for all of those genomes
    ani(
        clean_genomes.out
    )

}