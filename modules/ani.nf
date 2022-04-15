include {
    mash_dist;
    mash_sketch;
    mash_join;
    aggregate_distances;
} from './processes/ani'

include { clean_genomes } from './processes/align_genomes'

workflow ani {
    take:
        genomes_ch

    main:
        // Compute compressed genome representations (sketches) with mash
        mash_sketch(
            genomes_ch
        )

        // Combine all of the sketches
        mash_join(
            mash_sketch.out.toSortedList()
        )

        // Calculate the distance of each genome against the total
        mash_dist(
            mash_sketch.out.combine(mash_join.out)
        )

        // Join together those distances
        aggregate_distances(
            mash_dist.out.toSortedList()
        )

    emit:
    distances = aggregate_distances.out
}

workflow {

    // If the user specified the --genomes flag
    if ( params.genomes ) {

        // Then make a channel of the files specified there
        Channel
            .fromPath(params.genomes)
            .ifEmpty { error "Cannot find any files at ${params.genomes} -- use wildcard to specify all files in a folder" }
            .set { genomes_ch }

    // If the user did not specify the --genomes flag
    } else {

        // Then make a channel of the files in {project_folder}/genomes/
        Channel
            .fromPath("${params.project_folder}/genomes/*")
            .set { genomes_ch }

    }

    // Add any genomes found in the {project_folder}/downloaded_genomes/
    Channel
        .fromPath("${params.project_folder}/downloaded_genomes/*")
        .set { downloaded_genomes_ch }

    // Clean all of those genomes
    clean_genomes(
        genomes_ch
            .mix(downloaded_genomes_ch)
    )

    // Calculate the ANI for all of those genomes
    ani(
        clean_genomes.out
    )

}