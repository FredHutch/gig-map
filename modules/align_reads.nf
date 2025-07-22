// Import processes from processes/align_reads.nf

include {
    centroids_length;
    count_reads;
    diamond;
    concat_aln;
    concat_readcounts;
    diamond_logs;
    filter_aln;
    famli;
    famli_logs;
    gather
} from './processes/align_reads'

include { concat_csv as concat_diamond_logs } from './processes/general' addParams(output_csv: "diamond_logs.csv", output_subfolder: "logs", publish: true)

// Reuse the DIAMOND index creation process
include { makedb_diamond } from './processes/align_genomes'

workflow align_reads {
    take:
    // Single file containing the centroids in FASTA format
    centroids_faa
    // Channel of FASTQ files to align against those centroids
    fastq_ch

    main:

    // Get the length of each gene in the centroids
    // This is used to compute the RPKM later on
    centroids_length(
        centroids_faa
    )

    // Make a DIAMOND database for the centroids
    makedb_diamond(centroids_faa)

    // Count the number of reads provided for each specimen
    count_reads(
        fastq_ch
    )

    // Join the read counts from each sample which are split across lanes
    count_reads.out
        .groupTuple()
        .branch {
            i -> 
            multiple: i[1].size() > 1
            single: true
        }
        .set { count_reads_grouped }

    concat_readcounts(
        count_reads_grouped.multiple
    )

    // Align those reads with DIAMOND against the centroids
    diamond(
        makedb_diamond.out,
        fastq_ch
    )

    // Summarize the time of execution for each alignment
    diamond_logs(
        diamond.out.log
    )

    // Join the stats across all files
    concat_diamond_logs(
        diamond_logs.out.toSortedList()
    )

    // Join the alignments from each sample which are split across lanes
    // The items which are grouped together are the ones which have the same sample name
    // If there is a sample with just one lane, it will be left alone
    diamond.out
        .aln
        .groupTuple()
        .branch {
            i -> 
            multiple: i[1].size() > 1
            single: true
        }
        .set { aln_grouped }

    concat_aln(
        aln_grouped.multiple
    )

    // Filter out any samples which have an insufficient number of alignments
    filter_aln(
        aln_grouped
            .single
            .map {
                it -> [it[0], it[1][0]]
            }
            .mix(concat_aln.out)
    )

    // Filter the alignments and resolve multi-mapping reads with FAMLI
    famli(filter_aln.out)

    // Summarize the amount of time used for FAMLI
    famli_logs(
        famli.out.log.toSortedList()
    )

    // Gather all of the alignment information into a single CSV
    gather(
        famli.out.json.toSortedList(),
        count_reads_grouped
            .single
            .map {
                it -> it[1][0]
            }
            .mix(concat_readcounts.out)
            .toSortedList()
    )

    emit:
    csv = gather.out
    centroids_length = centroids_length.out

}