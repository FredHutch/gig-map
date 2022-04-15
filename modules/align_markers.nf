
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
    marker_genes

    main:
    
    // Build a BLAST database with the marker sequences
    makedb_blast(
        marker_genes
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
        "markers"
    )

    // Run the MSA
    combine_markers(
        reorganize_fastas.out.flatten()
    )

    // Generate phylogenetic trees for each marker gene
    raxml(
        combine_markers.out.msa
    )

    emit:
        distmat = combine_markers.out.distmat
        msa = combine_markers.out.msa
}