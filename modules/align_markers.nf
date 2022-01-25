
// Import the processes to run in this workflow
include {
    makedb_blast;
    align_blast;
    filter_alignments;
    extract_genes;
    reorganize_fastas;
} from './processes/align_genomes'

include { 
    raxml;
    combine_markers;
} from './processes/align_markers'

workflow align_markers {
    take:
    genomes_ch
    upstream_markers_ch

    main:

    // If the user specified the --markers flag/param
    if ( params.markers ) {

        // Read markers from that wildcard/path
        markers_fp = params.markers
    
    // If the user did not specify the --markers flag/param
    } else {

        // Read any files in the markers/ folder in the project folder
        markers_fp = "${params.project_folder}/markers/"
    }

    // Combine any files from those folders with the markers
    // selected by alignment to the genomes
    Channel
        .fromPath(markers_fp)
        .mix(upstream_markers_ch)
        .ifEmpty { error "No marker genes found" }
        .set { markers_ch }
    
    // Build a BLAST database with the marker sequences
    makedb_blast(
        markers_ch
            .toSortedList()
    )

    // Align the genomes against those markers
    align_blast(
        makedb_blast.out,
        genomes_ch
    )

    // Filter any overlapping alignments
    filter_alignments(
        align_blast.out
    )

    // Extract the aligned regions for each marker,
    // combining the markers provided by the user with the 
    // markers which were discovered from the input data
    extract_genes(
        filter_alignments.out.join(
            genomes_ch.map({
                it -> [it.name, it]
            })
        )
    )

    // The output of extract_markers has one file per genome, with the
    // FASTA headers indicating the marker of origin

    // Next we will reformat the markers to have one file per marker,
    // with the FASTA headers indicating the genome of origin
    reorganize_fastas(
        extract_genes.out.toSortedList(),
        "markers/fasta/"
    )

    // Run the MSA
    combine_markers(
        reorganize_fastas.out.flatten()
    )

    // Generate phylogenetic trees for each marker gene
    raxml(
        combine_markers.out.msa
    )
}