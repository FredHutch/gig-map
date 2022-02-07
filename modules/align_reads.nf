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

    main:

    // Decide which input folder to use
    single_end_folder = ( params.single_end_fastq ) ? params.single_end_fastq : "${params.project_folder}/metagenomes/single_end"
    paired_end_folder = ( params.paired_end_fastq ) ? params.paired_end_fastq : "${params.project_folder}/metagenomes/paired_end"

    // Make a channel with any single-end reads
    Channel
        .fromPath(
            "${single_end_folder}/*${params.single_end_suffix}"
        )
        .map {
            it -> [
                // Remove the suffix from the file name
                it.name.substring(0, it.name.length() - params.single_end_suffix.length()),
                // Include the file in a tuple with the modified name
                // Note that this file is being placed in a list, to match
                // the pattern of the paired-end reads below
                [it]
            ]
        }
        .set { single_end_ch }

    // Making a channel with paired-end reads is a little easier
    Channel
        .fromFilePairs("${paired_end_folder}/*${params.paired_end_suffix}")
        .set { paired_end_ch }
    
    // In case there is any sample which has both single-end and paired-end reads,
    // Combine both of those channels
    single_end_ch
        // First flatten in to a list of tuples with [name, fastq]
        .transpose()
        // Add in the paired-end data
        .mix(
            // After similarly flattening the paired-end reads
            paired_end_ch
                .transpose()
        )
        // Finally, group by the sample name again
        .groupTuple()
        // Assign to fastq_ch
        .set { fastq_ch }

    // Only run the alignment if there are any files to align
    if ( fastq_ch.toSortedList().length() > 0 ){

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

        csv = gather.out

    } else {
        
        csv = Channel.empty()

    }

    emit:
    csv = csv

}