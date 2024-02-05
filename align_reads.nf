#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { align_reads } from './modules/align_reads'
include { find_reads } from './modules/find_reads'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Align genes against a set of WGS sequence reads
        
        Aligns a deduplicated collection of genes against a collection of whole-genome shotgun
        short-read sequencing data, producing a summary of the number of reads which align to each gene

        Parameters:

        --genes             Path of single deduplicated amino acid FASTA file to be used for alignment
        --reads             Folder containing the set of FASTQ files to align against those genes
        --reads_suffix      File ending for all reads (default: ${params.reads_suffix})
        --paired            Set to true if input folder contains paired-end FASTQ files (default: false)
        --samplesheet       Optional file listing the FASTQ inputs to process (headers: sample,R1,R2)
        --output            Folder where output files will be written

        --min_score_reads   Minimum alignment score (default: ${params.min_score_reads})
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
    helpers.require_param(params.genes, "genes")

    // If there is no samplesheet
    if ( "${params.samplesheet}" == "false" ){
        // Requre the reads parameter
        helpers.require_param(params.reads, "reads")
    }

    // Get the reads either from samplesheet or the input folder
    find_reads()

    // Align those reads against the centroids
    align_reads(
        file("${params.genes}"),
        find_reads.out
    )

}