#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflow
include { render; serialize } from './modules/render'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Render genes-in-genomes map
        
        Generate a visual display summarizing the gene content of a set of genomes

        Parameters:

        --output            Folder where output files will be written
        --genome_aln        Genome alignments (CSV)
        --genome_annot      Genome annotations, keyed on the column 'genome_id' (CSV)
        --gene_annot        Gene annotations, keyed on the column 'gene_id' (CSV)
        --gene_order        TXT file containing the ordered list of genes
        --genome_distmat    Distance matrix
        --render_options    Configuration options for the gig-map display (JSON)

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output,         "output")
    helpers.require_param(params.genome_aln,     "genome_aln")
    helpers.require_param(params.genome_annot,   "genome_annot")
    helpers.require_param(params.gene_annot,     "gene_annot")
    helpers.require_param(params.gene_order,     "gene_order")
    helpers.require_param(params.genome_distmat, "genome_distmat")
    helpers.require_param(params.render_options, "render_options")

    // The input files all must exist
    genome_aln          = Channel.fromPath(params.genome_aln, checkIfExists: true)
    genome_annot        = Channel.fromPath(params.genome_annot, checkIfExists: true)
    gene_annot          = Channel.fromPath(params.gene_annot, checkIfExists: true)
    gene_order          = Channel.fromPath(params.gene_order, checkIfExists: true)
    genome_distmat      = Channel.fromPath(params.genome_distmat, checkIfExists: true)
    render_options      = Channel.fromPath(params.render_options, checkIfExists: true)

    // Run the collect sub-workflow
    render(
        genome_aln,
        genome_annot,
        genome_distmat,
        gene_order,
        gene_annot,
        render_options
    )

    // Save the files to disk
    serialize(
        genome_aln,
        genome_annot,
        genome_distmat,
        gene_order,
        gene_annot,
        render_options
    )

}