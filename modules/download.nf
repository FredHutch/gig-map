
// Workflow dedicated to download only

include { download_genomes } from './download_genomes'
include { download_genes } from './download_genes'

workflow download_wf {

    download_genomes()
    download_genes()

    emit:
    genomes = download_genomes.out.genomes
    genes = download_genes.out.genes
    
}
