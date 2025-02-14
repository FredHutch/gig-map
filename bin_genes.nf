#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { bin_genes_wf } from './modules/bin_genes'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Bin genes based on which genomes they are found within
        
        Processes the output of align_genomes.nf, grouping genes into bins and
        summarizing the genomes on the basis of those bins.

        Steps:
            1. Filter alignments by coverage and identity thresholds
            2. Filter genes by the number of genomes each is found in
            3. Group genes into bins
            4. Filter bins by minimum size (number of genes per bin)
            5. Group genomes by shared gene bin composition

        Parameters:

        --genome_aln        Alignments from align_genomes.nf (e.g. genomes.aln.csv.gz)
        --gene_annot        Annotations of each gene from deduplicate.nf (e.g. centroids.annot.csv.gz)
        --output            Folder where output files will be written

        --min_coverage      Minimum proportion of a gene which must align in order to retain each alignment
                            (default: ${params.min_coverage}, ranges 0-100)
        --min_identity      Minimum percent identity of the amino acid alignment required to retain each alignment
                            (default: ${params.min_identity}, ranges 0-100)

        --min_genomes_per_gene  Minimum number of genomes for a gene to be found in to be included (default: 1)

        --max_dist_genes    Maximum Jaccard distance threshold used to group genes into bins

        --min_bin_size      Minimum number of genes needed to retain a bin

        --max_dist_genomes  Maximum Euclidean distance threshold used to group genomes based on gene bin content


        Outputs:

        gene_bins.csv       Table listing genes and bins (including optional gene annotations)
        genome_groups.csv   Table listing the genome groups based on gene bin content
        bin_genes.log       Logfile reporting processing metrics
        *.html              Various summary figures
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genome_aln, "genome_aln")

    // Get the genome alignments
    genome_aln = file(params.genome_aln, checkIfExists: true)

    // Get the gene annotations
    gene_annot = file(params.gene_annot, checkIfExists: true)

    // Get the genome annotations
    genome_annot = file(params.genome_annot, checkIfExists: true)

    // Bin the genes
    bin_genes_wf(
        genome_aln,
        gene_annot,
        genome_annot
    )

}