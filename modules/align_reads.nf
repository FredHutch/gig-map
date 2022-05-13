// Import processes from processes/align_reads.nf

include {
    count_reads;
    diamond;
    famli;
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

    // Filter the alignments and resolve multi-mapping reads with FAMLI
    famli(diamond.out)

    // Gather all of the alignment information into a single CSV
    gather(
        famli.out.toSortedList(),
        count_reads.out.toSortedList()
    )

    emit:
    csv = gather.out

}