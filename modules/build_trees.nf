// Using DSL-2
nextflow.enable.dsl=2

include { 
    gene_msa;
    calc_distmat;
    merge_bins;
    raxml;
    filter_msas;
} from './processes/align_genes'

workflow build_trees {
    take:
        gene_fastas
        gene_bins
    main:
        // Align the gene FASTA files
        gene_msa(gene_fastas)

        // Calculate distance matrices from the MSAs
        calc_distmat(
            gene_msa.out.msa.ifEmpty { error "No MSAs were produced from the input gene FASTA files" }
        )

        // Combine the distance matrices for each bin
        merge_bins(
            gene_msa.out.msa.ifEmpty { error "No MSAs were produced from the input gene FASTA files" }.toSortedList(),
            calc_distmat.out.distmat.ifEmpty { error "No distance matrices were produced from the input gene FASTA files" }.toSortedList(),
            gene_bins
        )

        // Filter the MSAs based on the minimum thresholds for how much useful information is present,
        // both on a per-site and per-genome basis
        // This helps to avoid building trees from very sparse data
        filter_msas(
            merge_bins.out
            .msa
            .ifEmpty { error "No combined MSAs were produced from the merge_bins process" }
            .flatten()
        )

        // Build ML trees from the combined MSAs
        raxml(filter_msas.out)
}