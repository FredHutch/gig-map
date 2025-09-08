
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
    gene_msa;
    calc_distmat;
} from './processes/align_genes'

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
    gene_msa(
        reorganize_fastas.out.flatten()
    )

    // Generate the distance matrices
    calc_distmat(
        gene_msa.out.msa.ifEmpty { error "No MSAs were produced from the input gene FASTA files" }
    )

    // Generate phylogenetic trees for each marker gene
    raxml(
        gene_msa.out.msa.ifEmpty { error "No MSAs were produced from the input gene FASTA files" }
    )

    emit:
        distmat = calc_distmat.out.distmat
        msa = gene_msa.out.msa
}