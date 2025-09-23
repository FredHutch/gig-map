include {
    bin_summary;
    collect;
    split;
    plot
} from "./processes/bin_metagenomes"

include { regress } from "./test_reads"

workflow bin_metagenomes {
    take:
        read_alignments
        centroids_length
        gene_bins
        genome_groups
        group_profile
        metadata

    main:

        bin_summary(
            read_alignments,
            gene_bins,
            centroids_length
        )
        if (params.wide_metagenome_output){

            collect(
                read_alignments,
                gene_bins,
                genome_groups,
                group_profile,
                metadata
            )

            if ("${params.formula}" != "false"){

                regress(
                    collect.out.metadata,
                    collect.out.bin_counts
                )

                split(
                    regress.out
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
        }

}