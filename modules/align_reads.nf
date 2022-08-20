// Import processes from processes/align_reads.nf

include {
    count_reads;
    diamond;
    gather_logs;
    filter_aln;
    famli;
    famli_logs;
    gather
} from './processes/align_reads'

// Reuse the DIAMOND index creation process
include { makedb_diamond } from './processes/align_genomes'

workflow align_reads {
    take:
    // Single file containing the centroids in FASTA format
    centroids_faa
    // Channel of FASTQ files to align against those centroids
    fastq_ch

    main:

    // Make a DIAMOND database for the centroids
    makedb_diamond(centroids_faa)

    // Count the number of reads provided for each specimen
    count_reads(
        fastq_ch
    )

    // Align those reads with DIAMOND against the centroids
    diamond(
        makedb_diamond.out,
        fastq_ch
    )

    // Summarize the time of execution for each alignment
    gather_logs(
        diamond.out.log.toSortedList()
    )

    // Filter out any samples which have an insufficient number of alignments
    filter_aln(
        diamond.out.aln
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
        count_reads.out.toSortedList()
    )

    emit:
    csv = gather.out

}