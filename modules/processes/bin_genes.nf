process bin_genes {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
    path genome_aln
    path gene_annot

    output:
    path "*"

    script:
    template "bin_genes.py"
}