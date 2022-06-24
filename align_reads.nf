#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { align_reads } from './modules/align_reads'
include { join_read_pairs } from './modules/processes/align_reads'

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
    helpers.require_param(params.reads, "reads")
    helpers.require_param(params.genes, "genes")

    // If the input is paired-end
    if ( "${params.paired}" != "false" ) {

        // Get reads as pairs of files which differ only by containing '1' vs '2'
        Channel
                .fromFilePairs([
                    "${params.reads}**{1,2}*${params.reads_suffix}",
                    "${params.reads}**R{1,2}*${params.reads_suffix}",
                    "${params.reads}**R{1,2}_001${params.reads_suffix}",
                ])
                .ifEmpty { error "No reads found at ${params.reads}**{1,2}*${params.reads_suffix}"}
                .set { fastq_ch }

        join_read_pairs(fastq_ch)
        reads_ch = join_read_pairs.out

    } else {

        // Get reads as any file with the expected ending
        // and add the name (without the suffix) to match the output of join_read_pairs
        reads_ch = Channel
            .fromPath("${params.reads}**${params.reads_suffix}")
            .map {
                it -> [it.name.substring(0, it.name.length() - "${params.reads_suffix}".length()), it]
            }

    }

    // Align those reads against the centroids
    align_reads(
        file("${params.genes}"),
        reads_ch
    )

}