// Cluster genes by amino acid similarity
process cdhit {
    container "${params.container__cdhit}"
    label 'mem_medium'
    publishDir "${params.project_folder}/deduplicated_genes", mode: 'copy', overwrite: true
   
    input:
    path "input.genes.*.fasta.gz"
    
    output:
    path "centroids.faa.gz", emit: fasta
    path "centroids.membership.csv.gz"
    
    script:
    template "cdhit.sh"
}


// Generate a simple annotation path for each centroid
process annotate_centroids {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.project_folder}", mode: 'copy', overwrite: true
   
    input:
    path "clustered.genes.fasta.gz"
    
    output:
    path "clustered.genes.csv.gz", emit: annot
    
    script:
    template "annotate_centroids.py"
}
