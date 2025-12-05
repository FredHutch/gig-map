include {
    bin_summary;
    plot_regress;
    wide_bin_abundance
} from "./processes/bin_metagenomes"

include { regress } from "./test_reads"

include { assign_metadata } from "./processes/assign_metadata"

workflow contrast_metagenomes {
    take:
        read_alignments
        centroids_length
        gene_bins
        metadata

    main:

        bin_summary(
            read_alignments,
            gene_bins,
            centroids_length
        )

        wide_bin_abundance(bin_summary.out.bin_summary)

        assign_metadata(metadata)

        regress(
            assign_metadata.out,
            wide_bin_abundance.out.fragments_per_million
        )

        plot_regress(regress.out.results)

}