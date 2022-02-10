
// Import the processes to run in this workflow
include {
    clean_genomes;
    filter_genes;
    makedb_blast;
    align_blast;
    align_diamond;
    filter_alignments;
    makedb_diamond;
    add_genome_name;
    concatenate_alignments;
    order_genes;
    select_markers;
    extract_genes;
    reorganize_fastas;
} from './processes/align_genomes'

workflow align_genomes {

    take:
    input_genomes_ch
    centroids_faa

    main:

    // If the --genomes param has been set by the user
    if ( params.genomes ) {

        // Get all of the files which are specified
        Channel
            .fromPath(params.genomes)
            .mix(input_genomes_ch)
            .set { genomes_ch }

    // Otherwise, if the --genomes param has not been set
    } else {

        // Read files from the genomes/ subfolder in the project folder
        Channel
            .fromPath("${params.project_folder}/genomes/*.gz")
            .mix(input_genomes_ch)
            .set { genomes_ch }

    }

    // Clean up the genome formatting
    clean_genomes(
        genomes_ch
            .map {
                it -> [it.name, it]
            }
    )

    // If the user has selected DIAMOND for alignment
    if (params.aligner == "diamond"){

        // Make a DIAMOND database for the input genes
        makedb_diamond(
            centroids_faa
        )

        // Align the query genes against the genomes
        align_diamond(
            makedb_diamond.out,
            clean_genomes.out
        )

        // Channel with all alignment results
        alignments_output = align_diamond.out
    }

    // If the user has selected BLAST for alignment
    if (params.aligner == "blast"){

        // Make a BLAST database for the input genes
        makedb_blast(
            centroids_faa
        )

        // Align the query genes against the genomes
        align_blast(
            makedb_blast.out,
            clean_genomes.out
        )

        // Channel with all alignment results
        alignments_output = align_blast.out
    }

    // Filter any overlapping alignments
    filter_alignments(
        align_diamond.out
    )

    // Add the name of the query genome to the alignments file
    add_genome_name(
        filter_alignments.out
    )

    // Concatenate the results
    concatenate_alignments(
        add_genome_name.out.toSortedList()
    )

    // Select a set of marker genes from the provided alignments
    select_markers(
        concatenate_alignments.out,
        centroids_faa
    )

    // Order the genes based on the genomes they align to
    order_genes(
        concatenate_alignments.out
    )

    // Extract the aligned regions for all genes
    extract_genes(
        alignments_output.join(
            clean_genomes.out.map({
                it -> [it.name, it]
            })
        )
    )

    // Next we will reformat the genes to have one file per gene,
    // with the FASTA headers indicating the genome of origin
    reorganize_fastas(
        extract_genes
            .out
            .ifEmpty { error "No genes found" }
            .toSortedList(),
        "genome_alignments/genes/"
    )

    emit:
    markers = select_markers.out
    clean_genomes = clean_genomes.out
    concat_alignments = concatenate_alignments.out
    gene_order = order_genes.out

}