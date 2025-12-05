include {
    bin_summary;
    collect;
    split;
    plot_metagenomes;
    plot_regress;
    wide_bin_abundance
} from "./processes/bin_metagenomes"

include { regress } from "./test_reads"

include { assign_metadata } from "./processes/assign_metadata"

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

            wide_bin_abundance(bin_summary.out.bin_summary)

            collect(
                read_alignments,
                gene_bins,
                genome_groups,
                group_profile,
                metadata
            )

            if ("${params.formula}" != "false"){

                assign_metadata(collect.out.metadata)

                regress(
                    assign_metadata.out,
                    wide_bin_abundance.out.fragments_per_million
                )

                plot_regress(regress.out.results)

                split(
                    regress.out.results
                )

                plot_metagenomes(
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