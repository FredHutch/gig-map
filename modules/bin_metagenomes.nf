include {
    collect;
    split;
    plot
} from "./processes/bin_metagenomes"

include { corncob } from "./test_reads"

workflow bin_metagenomes {
    take:
        read_alignments
        gene_bins
        genome_groups
        group_profile
        metadata

    main:

        collect(
            read_alignments,
            gene_bins,
            genome_groups,
            group_profile,
            metadata
        )

        corncob(
            collect.out.metadata,
            collect.out.bin_counts
        )

        split(
            corncob.out
        )

        plot(
            split
                .out
                .flatten()
                .map {
                    it -> [it.name.replaceAll(".results.csv", ""), it]
                }
                .combine(collect.out.bins_h5ad)
        )

}