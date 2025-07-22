#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { bin_metagenomes } from './modules/bin_metagenomes'
include { align_reads } from './modules/align_reads'
include { find_reads } from './modules/find_reads'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Align a set of metagenomic sequencing data to a binned pangenome in
        which the genes have been grouped by similar occurrence across genomes
        (using the build_pangenome.nf workflow).

        Those binned gene groups will be used to summarize a batch of metagenomic
        data, and the similarity of gene abundance across metagenomes will be
        tested for concordence with the similarity across genomes.

        By applying the results of gene binning to the metagenome alignment
        results using the same gene catalog, the diversity of a collection of
        organisms can be collapsed into units of those co-occurrent genetic
        elements.

        Parameters:

        --genes             Path of single deduplicated amino acid FASTA file to be used for alignment
        --reads             Folder containing the set of FASTQ files to align against those genes
        --reads_suffix      File ending for all reads (default: ${params.reads_suffix})
        --paired            Set to true if input folder contains paired-end FASTQ files (default: false)
        --samplesheet       Optional file listing the FASTQ inputs to process (headers: sample,R1,R2)

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

        --gene_bins         Grouping of genes into bins (e.g. gene_bins.csv)
        --group_profile     Gene content of genome groups (e.g. group_profile.csv)
        --genome_groups     Grouping of genomes by gene content (e.g. genome_groups.csv)
        --output            Folder where output files will be written

        --metadata          Optional: Metadata table used to compare samples (CSV)
        --formula           Optional: Column(s) from metadata table used for comparison

        --min_n_reads       Exclude any samples with fewer than this number of reads aligned (default: 0)
        --min_n_genes       Exclude any samples with fewer than this number of genes detected (default: 0)

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.genes, "genes")
    helpers.require_param(params.gene_bins, "gene_bins")
    helpers.require_param(params.group_profile, "group_profile")
    helpers.require_param(params.output, "output")

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

    // Get the gene bins
    gene_bins = file(params.gene_bins, checkIfExists: true)

    // Get the genome groups
    genome_groups = file(params.genome_groups, checkIfExists: true)

    // Get the genome groups
    group_profile = file(params.group_profile, checkIfExists: true)

    // Get the metadata table
    metadata = file(params.metadata, checkIfExists: true)

    // Bin the metagenomes
    bin_metagenomes(
        align_reads.out.csv,
        align_reads.out.centroids_length,
        gene_bins,
        genome_groups,
        group_profile,
        metadata
    )

}