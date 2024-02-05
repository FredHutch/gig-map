include { join_read_pairs } from './processes/align_reads'

workflow find_reads {
    main:
    // If a samplesheet was provided
    if ( "${params.samplesheet}" != "false" ){
        Channel
            .from(
                file(
                    "${params.samplesheet}",
                    checkIfExists: true
                )
            )
            .splitCsv(
                header: true,
                strip: true
            )
            .map {
                it -> [
                    it.get('sample', it['"sample"']),
                    [
                        file(it.get('R1', it.get('fastq_1', it['"fastq_1"']).replaceAll(/^\"|\"$/, "")), checkIfExists: true),
                        file(it.get('R2', it.get('fastq_2', it['"fastq_2"']).replaceAll(/^\"|\"$/, "")), checkIfExists: true)
                    ]
                ]
            }
            .set { fastq_ch }

            join_read_pairs(fastq_ch)
            reads_ch = join_read_pairs.out

    } else {

        // Remove any trailing slash from the reads folder
        reads_folder = params.reads.replaceAll('/$', '')

        // If the input is paired-end
        if ( "${params.paired}" != "false" ) {

            // Get reads as pairs of files which differ only by containing '1' vs '2'
            Channel
                    .fromFilePairs(
                        "${reads_folder}/**${params.read_pairing_pattern}*${params.reads_suffix}"
                    )
                    .ifEmpty { error "No reads found at ${reads_folder}/**${params.read_pairing_pattern}*${params.reads_suffix}, consider modifying --reads, --read_pairing_pattern, or --reads_suffix"}
                    .set { fastq_ch }

            join_read_pairs(fastq_ch)
            reads_ch = join_read_pairs.out

        } else {

            // Get reads as any file with the expected ending
            // and add the name (without the suffix) to match the output of join_read_pairs
            reads_ch = Channel
                .fromPath("${reads_folder}/**${params.reads_suffix}")
                .map {
                    it -> [it.name.substring(0, it.name.length() - ("${params.reads_suffix}".length() + 1)), it]
                }

        }
    }

    emit:
    reads_ch
}