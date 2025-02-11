process filter_genomes_ids {
    container "${params.container__pandas}"
    label 'io_limited'

    input:
        path "distances.csv.gz"

    output:
        path "filtered_genomes_ids.txt"

    script:
    template "filter_genomes_ids.py"

}

workflow filter_genomes {
    take:
        distances
        genomes_ch

    main:
        // Filter genomes based on ANI
        // This process will emit a list of IDs
        filter_genomes_ids(distances)

        // Split the output text file from filter_genome_ids
        // into a list, splitting on the \n character
        filter_genomes_ids
            .out
            .splitCsv(sep: '\n')
            .set { filter_genomes_ids_ch }

        // Filter the genomes_ch to only include
        // those files which are in the list of IDs
        genomes_ch
            .map { [it.name, it] }
            .join(filter_genomes_ids_ch)
            .map { it[1] }
            .set { filtered_genomes }

    emit:
        filtered_genomes = filtered_genomes

}