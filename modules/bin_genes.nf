include { bin_genes; plot_bins } from './processes/bin_genes'

workflow bin_genes_wf {
    take:
        genome_aln
        gene_annot
        genome_annot

    main:
        bin_genes(
            genome_aln,
            gene_annot,
            genome_annot
        )

        plot_bins(
            genome_aln,
            bin_genes.out
        )

}
