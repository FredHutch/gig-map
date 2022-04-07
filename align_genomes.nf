#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflow
include { align_genomes } from './modules/align_genomes'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Align genes against a set of genomes
        
        Aligns a deduplicated collection of genes against a collection of genomes, producing
        a summary of the location of each gene alignment across each genome

        Parameters:

        --genes             Path of single deduplicated amino acid FASTA file to be used for alignment
        --genomes           Folder containing the set of genome nucleotide FASTAs to align against
        --output            Folder where output files will be written

        --min_coverage      Minimum proportion of a gene which must align in order to retain the alignment
                            (default: ${params.min_coverage}, ranges 0-100)
        --min_identity      Minimum percent identity of the amino acid alignment required to retain the alignment
                            (default: ${params.min_identity}, ranges 0-100)
        --max_evalue        Maximum E-value threshold used to filter all alignments
                            (default: ${params.max_evalue})
        --aligner           Algorithm used for alignment (default: ${params.aligner}, options: diamond, blast)
        --query_gencode     Genetic code used for conceptual translation of genome sequences
                            (default: ${params.query_gencode})
        --max_overlap       Any alignment which overlaps a higher-scoring alignment by more than this
                            amount will be filtered out (default: ${params.max_overlap}, range: 0-100)
        --aln_fmt           Column headings used for alignment outputs (see DIAMOND documentation for details)
                            (default: ${params.aln_fmt})
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genomes, "genomes")
    helpers.require_param(params.genes, "genes")

    // Remove any trailing slash from the genome folder
    genome_folder = params.genomes.replaceAll('/$', '')

    // Get all of the genomes
    Channel
        .fromPath("${genome_folder}/*")
        .ifEmpty { error "Cannot find any files at ${genome_folder}/*" }
        .set { genomes_ch }

    // Get the gene FASTA
    genes_faa = file(params.genes, checkIfExists: true)

    // Run the genome alignment sub-workflow
    align_genomes(
        genomes_ch,
        genes_faa
    )

}