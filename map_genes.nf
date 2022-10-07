#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import the process to run
include {
    shard;
    map_genes_blast;
    map_genes_diamond;
    filter;
    join
} from './modules/processes/map_genes'

// Standalone entrypoint
workflow {

    helpers.help_message(
        """
        Compare two collections of genes, finding the best match for each of the
        genes in the query set from the collection of genes in the reference set.
        In the results which are provided, each query gene will have no more than
        one alignment reported. However, each reference gene may have alignments
        reported to more than one query gene.

        The output file will be written in CSV format as gene_mapping.csv.gz in
        the output directory

        Parameters:

        --queries           Query gene collection (FASTA)
        --references        Reference gene collection (FASTA)
        --output            Output directory

        --min_coverage      Minimum proportion of a gene which must align in order to retain the alignment
                            (default: ${params.min_coverage}, ranges 0-100)
        --min_identity      Minimum percent identity of the amino acid alignment required to retain the alignment
                            (default: ${params.min_identity}, ranges 0-100)
        --max_evalue        Maximum E-value threshold used to filter all alignments
                            (default: ${params.max_evalue})
        --map_batchsize     Number of genes to align in a batch
                            (default: ${params.map_batchsize})
        --aligner           Algorithm used for alignment (default: ${params.aligner}, options: diamond, blast)
        --aln_fmt           Column headings used for alignment outputs (see DIAMOND documentation for details)
                            (default: ${params.aln_fmt})

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.queries, "queries")
    helpers.require_param(params.references, "references")

    // Shard the genes
    shard(
        file(params.queries, checkIfExists: true, glob: false)
    )

    // Run the alignment
    if ( "${params.aligner}" == "blast" ){
        map_genes_blast(
            shard.out.flatten(),
            file(params.references, checkIfExists: true, glob: false)
        )
        unfiltered_aln = map_genes_blast.out
    }else{
        if ( "${params.aligner}" == "diamond" ){
            map_genes_diamond(
                shard.out.flatten(),
                file(params.references, checkIfExists: true, glob: false)
            )
            unfiltered_aln = map_genes_diamond.out
        }else{
            error "Parameter 'aligner' must be diamond or blast, not ${params.aligner}"
        }
    }

    // Filter down to just the top alignment per query
    filter(unfiltered_aln)

    // Join the shards
    join(filter.out.toSortedList())

}