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
    join_pdist;
    collect_pdist as batch_pdist;
    collect_pdist;
} from './modules/processes/map_genes'

// Standalone entrypoint
workflow {

    helpers.help_message(
        """
        Find the amino acid similarities for all pairwise comparisons of genes
        in a collection. A minimum threshold will be set for the coverage of
        the alignment over the length of the genes.

        The output file will be written in JSON format as gene_pdist.json.gz in
        the output directory. The format of that file will be:

        {
            "gene1": {
                "gene2": 95,
                "gene3": 36
            },
            "gene2": {
                "gene3": 42
            }
        }

        Parameters:

        --genes             Query gene collection (FASTA)
        --output            Output directory

        --min_coverage      Minimum proportion of a gene which must align in order to retain the alignment
                            (default: ${params.min_coverage}, ranges 0-100)
        --max_evalue        Maximum E-value threshold used to filter all alignments
                            (default: ${params.max_evalue})
        --map_batchsize     Number of genes to align in a batch
                            (default: ${params.map_batchsize})
        --aligner           Algorithm used for alignment (default: ${params.aligner}, options: diamond, blast)

        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")
    helpers.require_param(params.genes, "genes")

    // Shard the genes
    shard(
        file(params.genes, checkIfExists: true, glob: false)
    )

    // Make all pairwise combinations
    shard
        .out
        .flatten()
        .combine(shard.out.flatten())
        .filter { it[0].name <= it[1].name }
        .set { permuted }

    // Run the alignment
    if ( "${params.aligner}" == "blast" ){
        map_genes_blast(permuted)
        aln = map_genes_blast.out
    }else{
        if ( "${params.aligner}" == "diamond" ){
            map_genes_diamond(permuted)
            aln = map_genes_diamond.out
        }else{
            error "Parameter 'aligner' must be diamond or blast, not ${params.aligner}"
        }
    }

    // Join the shards
    join_pdist(aln.collate(100))
    batch_pdist(join_pdist.out.collate(100))
    collect_pdist(batch_pdist.out.collect())
}