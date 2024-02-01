#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import helpers
GroovyShell shell = new GroovyShell()
def helpers = shell.parse(new File("${workflow.projectDir}/helpers.gvy"))

// Import sub-workflows
include { align_genomes } from './modules/align_genomes' addParams(output: "${params.output}/align")
include { download_genes } from './modules/download_genes'
include { download_genomes } from './modules/download_genomes'
include { deduplicate } from './modules/deduplicate' addParams(output: "${params.output}/gene_catalog")
include { bin_genes } from './modules/processes/bin_genes' addParams(output: "${params.output}/bin_pangenome")

// Standalone entrypoint
workflow {

    // Show help message if the user specifies the --help flag at runtime
    helpers.help_message(
        """
        Build a pangenome dataset from a collection of genes and genomes

        1. Genomes and genes will be downloaded from NCBI if needed
        2. Genes will be deduplicated to yield a nonredundant gene catalog
        3. Genomes will be aligned against the gene catalog
        4. Genes and genomes will be binned by co-occurrence

        Parameters:

        --genome_tsv        (optional) Table of genomes to download from NCBI (e.g. https://www.ncbi.nlm.nih.gov/datasets/genome/)
        --genome_csv        (optional) Table of genomes to download from legacy NCBI (e.g. https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/)
        --genes             (optional) Folder containing the set of gene amino acid FASTAs to deduplicate
        --genomes           (optional) Folder containing the set of genome nucleotide FASTAs to align against
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

        --min_genomes_per_gene  Minimum number of genomes for a gene to be found in to be included (default: 1)
        --max_dist_genes    Maximum Jaccard distance threshold used to group genes into bins
        --min_bin_size      Minimum number of genes needed to retain a bin
        --max_dist_genomes  Maximum Euclidean distance threshold used to group genomes based on gene bin content
        """,
        params.help
    )

    // Make sure that the required parameters were provided
    helpers.require_param(params.output, "output")

    //////////////////////////////////////
    // MANUALLY INPUT GENES AND GENOMES //
    //////////////////////////////////////

    // Remove any trailing slash from the genome and gene folder
    if(params.genomes){params.genomes = params.genomes.replaceAll('/$', '')}
    if(params.genes){params.genes = params.genes.replaceAll('/$', '')}

    // Get all of the genomes
    Channel
        .fromPath("${params.genomes}/**.f*")
        .filter { it.exists() }
        .set { genomes_ch }

    // Get all of the genes
    Channel
        .fromPath("${params.genes}/**.f*")
        .filter { it.exists() }
        .set { genes_ch }

    ////////////////////////////////
    // DOWNLOAD GENOMES FROM NCBI //
    ////////////////////////////////

    // Parse the CSV path(s), if provided
    if(params.genome_csv.length() > 1){
        Channel
            .fromPath( params.genome_csv.split(',').toList() )
            .filter { it.exists() }
            .set { genome_manifest_csv }
    }else{
        Channel.empty().set { genome_manifest_csv }
    }

    // Parse the TSV path(s), if provided
    if(params.genome_tsv.length() > 1){
        Channel
            .fromPath( params.genome_tsv.split(',').toList() )
            .filter { it.exists() }
            .set { genome_manifest_tsv }
    }else{
        Channel.empty().set { genome_manifest_tsv }
    }

    // Download the genomes for each genome
    download_genomes(
        genome_manifest_csv,
        genome_manifest_tsv
    )

    //////////////////////////////
    // DOWNLOAD GENES FROM NCBI //
    //////////////////////////////

    // Download the genes for each genome
    download_genes(
        genome_manifest_csv,
        genome_manifest_tsv
    )

    //////////////////////////////
    // DEDUPLICATE GENE CATALOG //
    //////////////////////////////

    // Combine genes from NCBI and from the folder input
    deduplicate(genes_ch.mix(download_genes.out.genes))

    //////////////////////
    // ALIGN TO GENOMES //
    //////////////////////

    // Run the genome alignment sub-workflow
    align_genomes(
        genomes_ch.mix(download_genomes.out.genomes),
        deduplicate.out.fasta
    )

    ///////////////////////////
    // BIN GENES AND GENOMES //
    ///////////////////////////

    // Bin the genes
    bin_genes(
        align_genomes.out.concat_alignments,
        deduplicate.out.annot
    )

}