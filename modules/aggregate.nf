include {
    aggregate_results;
    cluster_genomes;
    cluster_genomes as cluster_genomes_by_marker;
    create_genome_manifest;
    create_gene_manifest;
} from './processes/aggregate'

workflow aggregate {

    take:
        concat_alignments
        ani_distances
        markers_distmat
        gene_order

    main:

        // Make genome groups based on the percent identity
        // in the region of the marker genes
        cluster_genomes(
            concat_alignments.combine(
                ani_distances
            ).combine(
                Channel.of(
                    params.ani_thresholds.split(/,/)
                )
            )
        )

        // Make genome groups based on the percent identity
        // in each of the marker genes
        cluster_genomes_by_marker(
            concat_alignments.combine(
                markers_distmat
            ).combine(
                Channel.of(
                    params.ani_thresholds.split(/,/)
                )
            )
        )

        // Make an empty annotation file for the genomes
        create_genome_manifest(
            concat_alignments
        )

        // Make an empty annotation file for the genes
        create_gene_manifest(
            concat_alignments
        )

        // Group together all results into a single HDF5 file object
        aggregate_results(
            concat_alignments,
            gene_order,
            ani_distances,
            cluster_genomes.out.toSortedList(),
            cluster_genomes_by_marker.out.toSortedList(),
            markers_distmat.toSortedList()
        )

    emit:
    genome_annot = create_genome_manifest.out
    gene_annot = create_gene_manifest.out

}