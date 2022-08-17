#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflow
include { collect } from './modules/collect'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Collect gene alignments
        
        Aggregate gene alignment information for rapid downstream visualization, including
        genome-genome ANI calculation, marker gene alignment, and compiling an .rdb object.

        Parameters:

        --output            Folder where output files will be written
        --genomes           Wildcard path indicating a set of genome nucleotide FASTAs
        --genome_aln        CSV file containing genome alignments
        --marker_genes      FASTA file containing marker gene sequences

        --query_gencode     Genetic code used for conceptual translation of genome sequences
                            (default: ${params.query_gencode})
        --max_evalue        Maximum E-value threshold used to filter marker alignments
                            (default: ${params.max_evalue})
        --aln_fmt           Column headings used for marker alignment outputs
                            (default: ${params.aln_fmt})
        --sketch_size       Number of minmers to use calculating genome-genome ANI distance
                            (default: ${params.sketch_size})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genomes, "genomes")
    helpers.require_param(params.marker_genes, "marker_genes")

    // Remove any trailing slash from the genome folder
    genome_folder = params.genomes.replaceAll('/$', '')

    // Get all of the genomes
    Channel
        .fromPath("${genome_folder}/*")
        .ifEmpty { error "Cannot find any files at ${genome_folder}/*" }
        .set { genomes_ch }

    // The genome alignment file must exist
    genome_aln = Channel.fromPath(params.genome_aln, checkIfExists: true)

    // The marker genes must also exist
    marker_genes = file(params.marker_genes, checkIfExists: true)

    // Run the collect sub-workflow
    collect(
        genomes_ch,
        genome_aln,
        marker_genes
    )

}