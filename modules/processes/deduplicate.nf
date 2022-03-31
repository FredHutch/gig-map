// Cluster genes by amino acid similarity
process cdhit {
    container "${params.container__cdhit}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
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
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path "clustered.genes.fasta.gz"
    
    output:
    path "centroids.annot.csv.gz", emit: annot
    
    script:
    template "annotate_centroids.py"
}
