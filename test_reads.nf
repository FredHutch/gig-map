#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { 
    parse_manifest;
    shard_genes;
    regress;
    join
} from './modules/test_reads'

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Test for differences in the relative abundance of genes between groups of specimens
        
        Uses the rigr:regress algorithm to test for significant differences in the proportional
        abundance of genes between different groups of samples.

        Parameters:

        --gene_abund        Read alignment data output by the `align_reads` tool
        --manifest          Table describing the experimental design, with `specimen` values matching the gene abundance table (CSV)
        --formula           Formula to test, using manifest column names as variables (e.g. GROUP, or TREATMENT + GROUP)
        --output            Folder where output files will be written

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.gene_abund, "gene_abund")
    helpers.require_param(params.manifest, "manifest")
    helpers.require_param(params.formula, "formula")
    helpers.require_param(params.output, "output")

    // Get the file with the gene abundances
    gene_abund = file(
        "${params.gene_abund}",
        checkIfExists: true,
        type: "file",
        glob: false
    )

    // Get the file with the manifest
    manifest = file(
        "${params.manifest}",
        checkIfExists: true,
        type: "file",
        glob: false
    )

    // Parse the manifest
    parse_manifest(
        manifest,
        gene_abund
    )

    // Shard the genes
    shard_genes(
        parse_manifest.out,
        gene_abund
    )

    // Test each gene
    regress(
        parse_manifest.out,
        shard_genes
            .out
            .flatten()
    )

    // Join the results
    join(
        regress.out.toSortedList()
    )

}